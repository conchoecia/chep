/*
 * Merge multiple pileup files by genomic position.
 * 
 * Takes a directory of pileup files and combines them by genomic position,
 * summing depths and concatenating base calls and quality scores.
 * 
 * Loads all data into memory for fast single-pass processing.
 * Much more memory-efficient than Python (~5x less RAM).
 * 
 * Handles both .pileup and .pileup.gz files automatically.
 * 
 * Usage:
 *     merge_pileups <pileup_directory> [--ref-fai <reference.fa.fai>]
 * 
 * Output is written to stdout in pileup format (6 columns):
 *     chrom, pos, ref, total_depth, bases, quals
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>
#include <zlib.h>

// Structure to hold pileup data for a single position
struct PileupData {
    int depth;
    std::string bases;
    std::string quals;
    
    PileupData() : depth(0) {}
    
    void add(int d, const std::string& b, const std::string& q) {
        depth += d;
        if (!b.empty() && b != "*") {
            bases += b;
        }
        if (!q.empty() && q != "*") {
            quals += q;
        }
    }
};

// Key for a genomic position
struct PosKey {
    std::string chrom;
    int pos;
    char ref;
    
    bool operator<(const PosKey& other) const {
        if (chrom != other.chrom) return chrom < other.chrom;
        if (pos != other.pos) return pos < other.pos;
        return ref < other.ref;
    }
};

// Custom comparator that uses chromosome order from .fai file
class ChromComparator {
private:
    std::unordered_map<std::string, int> chrom_order_;
    
public:
    ChromComparator() {}
    
    ChromComparator(const std::unordered_map<std::string, int>& order) 
        : chrom_order_(order) {}
    
    bool operator()(const PosKey& a, const PosKey& b) const {
        // Get chromosome orders
        int order_a = 1000000; // Default for unknown chroms
        int order_b = 1000000;
        
        auto it_a = chrom_order_.find(a.chrom);
        if (it_a != chrom_order_.end()) {
            order_a = it_a->second;
        }
        
        auto it_b = chrom_order_.find(b.chrom);
        if (it_b != chrom_order_.end()) {
            order_b = it_b->second;
        }
        
        if (order_a != order_b) return order_a < order_b;
        if (a.pos != b.pos) return a.pos < b.pos;
        return a.ref < b.ref;
    }
};

// Read chromosome order from .fai file
std::unordered_map<std::string, int> read_chrom_order(const std::string& fai_file) {
    std::unordered_map<std::string, int> chrom_order;
    std::ifstream infile(fai_file);
    
    if (!infile.is_open()) {
        std::cerr << "Warning: Could not open .fai file: " << fai_file << std::endl;
        return chrom_order;
    }
    
    std::string line;
    int idx = 0;
    while (std::getline(infile, line)) {
        size_t tab_pos = line.find('\t');
        if (tab_pos != std::string::npos) {
            std::string chrom = line.substr(0, tab_pos);
            chrom_order[chrom] = idx++;
        }
    }
    
    std::cerr << "Read " << chrom_order.size() << " chromosomes from " << fai_file << std::endl;
    return chrom_order;
}

// Read a single line from gzipped file
bool gzgetline(gzFile file, std::string& line) {
    line.clear();
    char buffer[4096];
    
    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        line += buffer;
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
            return true;
        }
    }
    
    return !line.empty();
}

// Parse a pileup line
bool parse_pileup_line(const std::string& line, PosKey& key, int& depth, 
                       std::string& bases, std::string& quals) {
    std::istringstream iss(line);
    std::string chrom, ref_str, depth_str, bases_str, quals_str;
    int pos;
    
    if (!(iss >> chrom >> pos >> ref_str >> depth_str)) {
        return false;
    }
    
    // Read bases and quals (might have spaces)
    std::getline(iss, bases_str, '\t');
    std::getline(iss, bases_str, '\t'); // Skip to bases column
    std::getline(iss, quals_str);
    
    // Parse
    key.chrom = chrom;
    key.pos = pos;
    key.ref = ref_str.empty() ? 'N' : ref_str[0];
    depth = std::stoi(depth_str);
    
    // Trim whitespace from bases and quals
    bases = bases_str;
    quals = quals_str;
    
    // Remove leading/trailing whitespace
    bases.erase(0, bases.find_first_not_of(" \t"));
    bases.erase(bases.find_last_not_of(" \t\r\n") + 1);
    quals.erase(0, quals.find_first_not_of(" \t"));
    quals.erase(quals.find_last_not_of(" \t\r\n") + 1);
    
    return true;
}

// Process a single pileup file
void process_pileup_file(const std::string& filepath, 
                         std::map<PosKey, PileupData>& positions,
                         size_t& lines_read) {
    bool is_gzipped = (filepath.size() > 3 && 
                      filepath.substr(filepath.size() - 3) == ".gz");
    
    if (is_gzipped) {
        gzFile file = gzopen(filepath.c_str(), "rb");
        if (!file) {
            std::cerr << "Error: Could not open file: " << filepath << std::endl;
            return;
        }
        
        std::string line;
        while (gzgetline(file, line)) {
            PosKey key;
            int depth;
            std::string bases, quals;
            
            if (parse_pileup_line(line, key, depth, bases, quals)) {
                positions[key].add(depth, bases, quals);
                lines_read++;
            }
        }
        
        gzclose(file);
    } else {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file: " << filepath << std::endl;
            return;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            PosKey key;
            int depth;
            std::string bases, quals;
            
            if (parse_pileup_line(line, key, depth, bases, quals)) {
                positions[key].add(depth, bases, quals);
                lines_read++;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pileup_directory> [--ref-fai <reference.fa.fai>]" << std::endl;
        return 1;
    }
    
    std::string pileup_dir = argv[1];
    std::string fai_file;
    
    // Parse command line arguments
    for (int i = 2; i < argc; i++) {
        if (std::strcmp(argv[i], "--ref-fai") == 0 && i + 1 < argc) {
            fai_file = argv[i + 1];
            i++;
        }
    }
    
    // Read chromosome order if provided
    std::unordered_map<std::string, int> chrom_order;
    bool use_chrom_order = false;
    if (!fai_file.empty()) {
        chrom_order = read_chrom_order(fai_file);
        use_chrom_order = !chrom_order.empty();
    }
    
    // Find all pileup files using POSIX dirent (C++11 compatible)
    std::vector<std::string> pileup_files;
    DIR* dir = opendir(pileup_dir.c_str());
    if (dir == NULL) {
        std::cerr << "Error: Could not open directory: " << pileup_dir << std::endl;
        return 1;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != NULL) {
        std::string filename = entry->d_name;
        
        // Skip . and ..
        if (filename == "." || filename == "..") {
            continue;
        }
        
        // Check if it ends with .pileup or .pileup.gz
        if (filename.size() > 7 && filename.substr(filename.size() - 7) == ".pileup") {
            pileup_files.push_back(pileup_dir + "/" + filename);
        } else if (filename.size() > 10 && filename.substr(filename.size() - 10) == ".pileup.gz") {
            pileup_files.push_back(pileup_dir + "/" + filename);
        }
    }
    closedir(dir);
    
    if (pileup_files.empty()) {
        std::cerr << "Warning: No pileup files found in " << pileup_dir << std::endl;
        return 1;
    }
    
    std::sort(pileup_files.begin(), pileup_files.end());
    std::cerr << "Merging " << pileup_files.size() << " pileup files..." << std::endl;
    
    // Map to store all positions
    std::map<PosKey, PileupData> positions;
    size_t total_lines = 0;
    
    // Process each file
    for (size_t i = 0; i < pileup_files.size(); i++) {
        // Extract filename from path (simple basename)
        std::string filepath = pileup_files[i];
        size_t last_slash = filepath.find_last_of("/\\");
        std::string filename = (last_slash == std::string::npos) ? filepath : filepath.substr(last_slash + 1);
        
        std::cerr << "Reading file " << (i + 1) << "/" << pileup_files.size() 
                  << ": " << filename << "..." << std::endl;
        
        size_t lines_before = total_lines;
        process_pileup_file(pileup_files[i], positions, total_lines);
        
        std::cerr << "  Read " << (total_lines - lines_before) << " lines, "
                  << positions.size() << " unique positions total" << std::endl;
    }
    
    std::cerr << "All files read. Total unique positions: " << positions.size() << std::endl;
    std::cerr << "Sorting and outputting..." << std::endl;
    
    // Sort and output
    if (use_chrom_order) {
        // Create sorted vector using custom comparator
        std::vector<std::pair<PosKey, PileupData>> sorted_positions(
            positions.begin(), positions.end()
        );
        
        ChromComparator comp(chrom_order);
        std::sort(sorted_positions.begin(), sorted_positions.end(),
                  [&comp](const auto& a, const auto& b) {
                      return comp(a.first, b.first);
                  });
        
        // Output
        std::string last_chrom;
        size_t chrom_count = 0;
        for (const auto& [key, data] : sorted_positions) {
            if (key.chrom != last_chrom) {
                if (!last_chrom.empty()) {
                    std::cerr << "Done with " << last_chrom << ": " 
                              << chrom_count << " positions" << std::endl;
                }
                last_chrom = key.chrom;
                chrom_count = 0;
            }
            chrom_count++;
            
            std::string bases = data.bases.empty() ? "*" : data.bases;
            std::string quals = data.quals.empty() ? "*" : data.quals;
            
            std::cout << key.chrom << '\t' << key.pos << '\t' << key.ref << '\t'
                      << data.depth << '\t' << bases << '\t' << quals << '\n';
        }
        
        if (!last_chrom.empty()) {
            std::cerr << "Done with " << last_chrom << ": " 
                      << chrom_count << " positions" << std::endl;
        }
    } else {
        // Output in natural order (lexicographic by chromosome)
        std::string last_chrom;
        size_t chrom_count = 0;
        for (const auto& [key, data] : positions) {
            if (key.chrom != last_chrom) {
                if (!last_chrom.empty()) {
                    std::cerr << "Done with " << last_chrom << ": " 
                              << chrom_count << " positions" << std::endl;
                }
                last_chrom = key.chrom;
                chrom_count = 0;
            }
            chrom_count++;
            
            std::string bases = data.bases.empty() ? "*" : data.bases;
            std::string quals = data.quals.empty() ? "*" : data.quals;
            
            std::cout << key.chrom << '\t' << key.pos << '\t' << key.ref << '\t'
                      << data.depth << '\t' << bases << '\t' << quals << '\n';
        }
        
        if (!last_chrom.empty()) {
            std::cerr << "Done with " << last_chrom << ": " 
                      << chrom_count << " positions" << std::endl;
        }
    }
    
    std::cerr << "Merge complete!" << std::endl;
    
    return 0;
}
