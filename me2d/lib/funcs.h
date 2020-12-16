#ifndef _FUNCS_H_
#define _FUNCS_H_

#include "defs.h"
#include <string>
#include <vector>

int check_quad();
int keyeq(const std::string &key, const std::string &str);
std::string string_lower(const std::string &str);
std::vector<std::string> split_string(const std::string &str, const std::string &delim);
std::string fmt_size(int64_t size);
double get_wtime();
std::string fmt_elapsed(const double t0, int add_parentheses=true);

#endif

