#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " sand.dat\n";
        return EXIT_FAILURE;
    }
    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Cannot open file: " << argv[1] << "\n";
        return EXIT_FAILURE;
    }

    std::string line;
    std::vector<int> prev_row;
    size_t width = 0;
    long line_no = 0;

    while (std::getline(infile, line)) {
        ++line_no;
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::vector<int> row;
        int h;
        // read heights until hitting a non-integer (g= token)
        while (ss >> h) {
            row.push_back(h);
        }
        if (row.empty()) {
            std::cerr << "Line " << line_no << ": no heights found\n";
            return EXIT_FAILURE;
        }
        if (width == 0) width = row.size();
        else if (row.size() != width) {
            std::cerr << "Line " << line_no << ": inconsistent row width (" << row.size() << " vs " << width << ")\n";
            return EXIT_FAILURE;
        }

        // Optionally skip parsing g= entirely

        // check neighbor invariant with periodic boundaries (except first line)
        if (!prev_row.empty()) {
            for (size_t i = 0; i < width; ++i) {
                int prev_h = prev_row[i];
                // skip check if height is 1
                if (prev_h == 1) continue;
                // periodic neighbor indices
                size_t left_idx  = (i == 0 ? width - 1 : i - 1);
                size_t right_idx = (i + 1 == width ? 0 : i + 1);
                int left = row[left_idx];
                int right = row[right_idx];
                if (left + right < prev_h) {
                    std::cerr << "Invariant violation at line " << line_no
                              << ", cell " << i
                              << ": left+right (" << (left + right)
                              << ") < prev_h (" << prev_h << ")\n";
                    return EXIT_FAILURE;
                }
            }
        }
        prev_row = std::move(row);
    }

    std::cout << "All invariants hold over " << line_no << " lines." << "\n";
    return EXIT_SUCCESS;
}