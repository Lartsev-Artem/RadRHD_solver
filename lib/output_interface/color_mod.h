/**
 * @file color_mod.h
 * @brief Управление цветом вывод в консоль
 * @version 0.1
 * @date 2024-02-09
 *
 */

#ifndef COLOR_MOD_H
#define COLOR_MOD_H

#include "ostream"

/**
 * @brief Пространство имён класса управляющего цветом вывода в консоль
 *
 */
namespace Color {
///\brief frontground color
enum eFG_Code {
  RESET = 0,
  FG_BLACK = 30,
  FG_RED = 31,
  FG_GREEN = 32,
  FG_YELLOW = 33,
  FG_BLUE = 34,
  FG_MAGENTA = 35,
  FG_CYAN = 36,
  FG_WHITE = 37,

  FG_DEFAULT = 39,
};

///\brief background color
enum eBG_Code {
  BG_BLACK = 40,
  BG_RED = 41,
  BG_GREEN = 42,
  BG_YELLOW = 43,
  BG_BLUE = 44,
  BG_MAGENTA = 45,
  BG_CYAN = 46,
  BG_WHITE = 47,

  BG_DEFAULT = 49
};

/**
 * @brief Класс изменяющий цвет вывода
 *
 */
class Modifier {
  eFG_Code fg_code; ///< цвет шрифта
  eBG_Code bg_code; ///< цвет фона
public:
  Modifier(eFG_Code fCode = FG_DEFAULT, eBG_Code gCode = BG_DEFAULT) : fg_code(fCode), bg_code(gCode) {}
  friend std::ostream &

  // изменяет весь вывод в консоль
  operator<<(std::ostream &os, const Modifier &mod) {
    return os << "\033[" << mod.bg_code << ";" << mod.fg_code << "m";
  }

  /// \brief Печатает конкретный текст указанным цветом
  void Printf(std::string str) {
    char last_symbol = ' ';
    if (str.back() == '\n') {
      str.back() = ' ';
      last_symbol = '\n';
    }

    printf("\033[%d;%dm%s\033[0m%c", bg_code, fg_code, str.c_str(), last_symbol);
  }
};
} // namespace Color
#endif //! COLOR_MOD_H