#include <string>
#include <fstream>
#include <vector>

/**
 * Reads all lines form a file and returns a vector of the lines.
 *
 * \param path Path to the file
 * \return A vector containing all the lines in the read file.
 */
std::vector<std::string> getlines(std::string path) {
  std::ifstream file(path);
  std::vector<std::string> strings;
  std::string line;
  while (std::getline(file, line)) {
    strings.push_back(line);
  }
  return strings;
}

