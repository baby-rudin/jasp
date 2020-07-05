#include "stropt.h"
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;

void    remove_comment(string &line)
{
    line += "!";
    line.erase(line.find_first_of("!"));
}

void    remove_space(string &line)
{
    line = " " + line;
    line.erase(0, line.find_first_not_of(" "));
    line.erase(line.find_last_not_of(" ") + 1);
}

void    clean_line(std::string &line)
{
    str_change(line, '\r', ' ');
    remove_comment(line);
    remove_space(line);
}

void    str_upper(std::string &str)
{ transform(str.begin(), str.end(), str.begin(), ::toupper); }

void    str_lower(std::string &str)
{ transform(str.begin(), str.end(), str.begin(), ::tolower); }

void    str_change(std::string &line, char pre, char aft)
{
    for (size_t i=0; i<line.size(); i++)
        if (line[i] == pre)
            line[i] = aft;
}

vector<string>  split_string(string line, string flag)
{
    vector<string>  ret;

    size_t pos = 0;
    std::string token;
    while ((pos = line.find(flag)) != std::string::npos) {
        token = line.substr(0, pos);
        ret.push_back(token);
        line.erase(0, pos + flag.length());
    }
    if (!line.empty())
        ret.push_back(line);

    return  ret;
}
