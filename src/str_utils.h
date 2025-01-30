#ifndef SURVEYOR_STRUTILS_H
#define SURVEYOR_STRUTILS_H

#include <vector>
#include <string>
#include <sstream>

//From Arafat Hasan: Answer to stack Overflow Question 14265581
std::vector<std::string> strsplit (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

#endif
