/**
 * @file parser_header.cpp
 * @brief Функция парсит файл с constexpr константами и формирует на его основе новый файл с вычисленными значениями f(x)
 * @date 2024-02-09
 *
 */
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "calculator.h"
#include "color_mod.h"

typedef double Type;      // базовый тип арифметических операций \warning значение (int=0) интерпретируется как внутрення ошибка и значение переменная не сохраняется
#define MATH_FUNC log     // используемая математическая функция
#define PREFIX_VAR "log_" // префикс у новых переменных

#define USE_ONLY_DOUBLE_TYPE // все типы приводить к double
//#define LOG_ENABLE         // включить лог

///\note при любой ситуации, не позволяющий сгенерировать конечный файл
#define ABORT abort() // exit(1)

#ifndef USE_ONLY_DOUBLE_TYPE
#include <map>
std::map<std::string, std::string>
    type_name_list;
#endif

static inline bool is_comments(std::string &str) {
  return (str.substr(0, 2) == "//");
}

static inline void skip_spaces(std::string &str) {
  int i = 0;
  while (isspace(str[i++])) {
  }
  str.erase(str.begin(), str.begin() + i - 1);
}

static inline bool skip_empty(std::ifstream &ifile, std::string &str) {
  bool is_end = false;
  while (str == "" && !is_end) {
    getline(ifile, str);
    is_end = ifile.eof();
  }
  return is_end;
}

static std::string get_file_name(const std::string &file) {
  const char *start = file.c_str();
  const char *end = start + file.length();

  while (*(end) != '.' && end != start)
    --end;

  const char *str = end;
  while (*str != '/' && *str != '\\' && (start != str))
    --str;
  str += (start != str); //пропускаем разделитель (если конец слова оставляем на месте)

#ifdef LOG_ENABLE
  std::cout << "Get input_file name: " << std::string(str, end) << "\n";
#endif
  return std::string(str, end);
}

static bool skip_section(std::ifstream &ifile, std::string &str) {
  static bool active_section = false;
  bool is_end = ifile.eof();

  if (str.find("#endif") != std::string::npos) {
    getline(ifile, str);
    return is_end;
  }

  if (str.find("#if 0") != std::string::npos) {
    active_section = false;
    do {
      getline(ifile, str);
    } while (str.find("#else") == std::string::npos);
    getline(ifile, str);
  }
  if (str.find("#if 1") != std::string::npos) {
    active_section = true;
    getline(ifile, str);
  }

  if (active_section && str.find("#else") != std::string::npos) {
    do {
      getline(ifile, str);
    } while (str.find("#endif") == std::string::npos);
    is_end = ifile.eof();
    getline(ifile, str);
    active_section = false;
  }

  return is_end;
}

static void parse_file(const std::string &file_input, Calculator<Type> &math) {
  std::ifstream ifile;
  ifile.open(file_input);
  if (!ifile.is_open()) {
    std::cout << "file: " << file_input << " wasn't opened\n";
    ABORT;
  } else {
    std::cout << "file: " << file_input << " has been opened!\n";
  }

  bool active_section = false;
  while (!ifile.eof()) {
    std::string str;
    std::getline(ifile, str);

    if (skip_empty(ifile, str)) {
      break;
    }

    skip_spaces(str);

    if (is_comments(str)) {
      continue;
    }

    if (skip_section(ifile, str)) {
      break;
    }

    const std::string cexpr = "constexpr";
    size_t pos_expr = str.find(cexpr);
    if (pos_expr == std::string::npos) {
      continue;
    }
    pos_expr += cexpr.length();

    std::string sub_str = str.substr(pos_expr, str.find(";") - pos_expr);

    skip_spaces(sub_str);

    int i = 0;
    while (sub_str[i++] != ' ') {
    }
    std::string type = sub_str.substr(0, i - 1);
    sub_str.erase(sub_str.begin(), sub_str.begin() + i);

#ifndef USE_ONLY_DOUBLE_TYPE
    skip_spaces(sub_str);

    i = 0;
    while (sub_str[i++] != ' ') {
    }
    std::string name = sub_str.substr(0, i - 1);

    type_name_list.insert(std::make_pair(name, type));
#endif

#ifdef LOG_ENABLE
    std::cout << "input string to calculator: " << sub_str << std::endl;
#endif

    Type n = math(sub_str.c_str());
    if (!n && math.was_error()) {
      std::cout << Color::Modifier(Color::FG_RED) << "Error while counting: " << math.error_message() << ": "
                << sub_str
                << Color::Modifier() << std::endl;
      ABORT;
    } else {
#ifdef LOG_ENABLE
      std::cout << "ans= " << n << std::endl;
#endif
    }
  }
  ifile.close();
}

static void write_file(const std::string &output_dir, const std::string &base_file, Calculator<Type> &math) {
  std::string name_f = PREFIX_VAR + base_file;
  std::string file_out = output_dir + name_f + ".h";
  std::ofstream ofile(file_out);
  if (!ofile.is_open()) {
    std::cout << "file: " << output_dir + name_f << " wasn't opened\n";
    ABORT;
  }

  std::transform(name_f.begin(), name_f.end(), name_f.begin(), [](unsigned char c) { return std::toupper(c); });

  ofile << "#ifndef " << name_f << "_H\n";
  ofile << "#define " << name_f << "_H\n";

  Calculator<Type>::itr_range rng;
  std::map<const char *, Type, Calculator<Type>::compcl>::iterator itr;

  rng = math.list_vars();
  for (itr = rng.begin; itr != rng.end; itr++) {
#ifndef USE_ONLY_DOUBLE_TYPE
    ofile << std::setprecision(16) << "constexpr " << type_name_list.find(itr->first)->second << " " << PREFIX_VAR << itr->first << " = " << (Type)MATH_FUNC(itr->second) << ";\n";
#else
    ofile << std::setprecision(16) << "constexpr double " << PREFIX_VAR << itr->first << " = " << MATH_FUNC(itr->second) << ";\n";
#endif
#ifdef LOG_ENABLE
    std::cout << itr->first << " = " << PREFIX_VAR << "(" << itr->second << ") -> " << (Type)MATH_FUNC(itr->second) << "\n";
#endif
  }

  ofile << "#endif //! " << name_f << "_H\n";
  ofile.close();

  std::cout << Color::Modifier(Color::FG_GREEN) << "file: " << file_out << " has been created!\n"
            << Color::Modifier();
}

int main(int argc, char **argv) {
  if (argc != 3) {
    printf("Error input format:\n");
    printf("input: input_file, output_dir\n");
    ABORT;
  }

  const std::string file_input = argv[1];
  const std::string output_dir = argv[2];

  // const std::string file_input = "test.h";
  // const std::string output_dir = "";

  Calculator<Type> math;
  parse_file(file_input, math);
  write_file(output_dir, get_file_name(file_input), math);

  return 0;
}
