
#ifndef STRING_UTILS_H
#define STRING_UTILS_H

template <typename str_type>
inline bool is_number(const str_type *str) {
  std::string s = *str;
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it))
    ++it;
  return !s.empty() && it == s.end();
}

#endif //! STRING_UTILS_H