#include "funcs.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <type_traits>
#include <limits>
#include <ios>
#include <cctype>
#include <ctime>
#ifdef _OPENMP
#include <omp.h>
#endif
using std::cout;
using std::endl;


int check_quad() {
  int digits10;
  if (std::is_same<quad_float, long double>::value) {
    digits10 = std::numeric_limits<long double>::digits10;
    if (digits10 < 33) {
      cout << "WARNING: check_quad: quadruple precision not enabled" << endl;
      cout << "         the program uses long double (" << digits10
        << " significant digits precision) instead." << endl;
    }
  }
  else { digits10 = 33; }
  return digits10;
}


int keyeq(const std::string &key, const std::string &str) {
  if (string_lower(key) == string_lower(str)) { return true; }
  return false;
}


std::string string_lower(const std::string &str) {
  int64_t i, n=str.length();
  std::string s;
  s.clear();
  for (i=0; i<n; i++) { s += (char)tolower(str[i]); }
  return s;
}

std::vector<std::string> split_string(const std::string &str, const std::string &delim) {
  int64_t i, n=str.length();
  std::vector<std::string> words;
  std::string word;
  word.clear();
  for (i=0; i<n; i++) {
    if (delim.find(str[i]) != std::string::npos) {
      if (!word.empty()) { words.push_back(word); }
      word.clear();
    }
    else { word += str[i]; }
  }
  if (!word.empty()) { words.push_back(word); }
  return words;
}

std::string fmt_size(int64_t size) {
  double dsize = (double)size / 1024.;
  std::ostringstream oss;
  std::string unit="KB";
  if (dsize >= 1000.) { dsize /= 1024.; unit = "MB"; }
  if (dsize >= 1000.) { dsize /= 1024.; unit = "GB"; }
  oss << std::fixed << std::setprecision(3) << dsize << " " << unit;
  return oss.str();
}


double get_wtime() {
  double t;
#ifdef _OPENMP
  t = omp_get_wtime();
#else
  t = (double)clock() / CLOCKS_PER_SEC;
#endif
  return t;
}

std::string fmt_elapsed(const double t0, int add_parentheses) {
  std::ostringstream oss;
  if (add_parentheses) oss << "(";
  oss << std::fixed << std::setprecision(1) << get_wtime() - t0 << " s";
  if (add_parentheses) oss << ")";
  return oss.str();
}



