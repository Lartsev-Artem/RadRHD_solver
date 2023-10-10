/**
 * @file utils.h
 * @brief Файл с настройками конфигурации сборки
 *
 * \details GENERATE_UTILS - макрос переключает цель сборки модуля.
    Если не определён, будут собрана статическая библиотека, доступная для
    вызова в других модулях.
    Если определён, то все имена функций переопределяются в main(argc argv)
    и генерируется набор исполняемых файлов

    \note для сборки требуется стандарт C++17

    \warning при сборке в библиотеку формат входных данных сохраняется,
    т.е. необходимо передать аргументы в формате стандартных аргументов командной строки.
    (первым элементом - имя программы)
 */
#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <string>

/*! \addtogroup utils Модуль утилит
    \brief Модуля включает в себя набор подпрограмм необходимых для
    пред- и пост- обработки данных
    @{
*/

//#define GENERATE_UTILS ///<генерировать исполняемые файлы из утилит

#ifdef GENERATE_UTILS
#define FUNC_NAME(_name) main
#else
#define FUNC_NAME(_name) utils::_name
#endif

/**
 * @brief Пространство имён модуля утилит
 *
 */
namespace utils {

template <typename T, typename T2>
void UtilsHelp(int argc, T str, T2 msg = "") {
  if (argc == 1) {
    return; //нет аргументов
  }

  std::string word = str[1];

  if (word.find("-h", 0, 2) == std::string::npos && word.find("--help", 0) == std::string::npos) {
    return;
  }

  setlocale(LC_CTYPE, "rus");
  std::cout << msg << '\n';
  exit(0);
}

#ifndef GENERATE_UTILS
/**
 * @brief Функция переводит нейтральный формат netgen в текстовый формат VTK
 *
 * @param[in] argc кол-во аргументов (требуется 3)
 * @param[in] argv массив char** содержит {файл_nethen[in] файл_vtk[out] разерность сетки [2/3]d  }
 * @return int ::e_type_completion
 */
int ReBuildNetgenToVTK(int argc, char **argv);

/**
 * @brief Функция создаёт копию сетки с добавлением нового скалярного поля
 *
 * @param argc[in] кол-во аргументов (требуется 3. 4 опционально)
 * @param argv[in] массив char** содержит{ файл_vtk[in] файл_данные.bin файл_vtk[out] [имя поля в структуре сетки]}
 * @return int ::e_type_completion
 */
int SetScalarDataVtkFromBinFile(int argc, char **argv);

/**
 * @brief Функция переводит текстовый файл с данными во внутренний формат
 *
 * @param argc[in] кол-во аргументов (требуется 2)
 * @param argv[in] массив char** содержит {файл данных[in] файл[out]}
 * @return int ::e_type_completion
 */
int DataToSimpleTxt(int argc, char **argv);

/**
 * @brief Функция добавляет нового скалярного поля в существующую сетку vtk
 *
 * @param argc[in] кол-во аргументов (требуется 2. 3 опционально)
 * @param argv[in] массив char** содержит { файл_vtk[in] файл_данные.bin [имя поля в структуре сетки]}
 * @return int ::e_type_completion
 */
int AddScalarDataVtkFromBinFile(int argc, char **argv);

/**
 * @brief Функция на основе бинарных файлов геометрии формирует два файла с геометрией в формате ::face_t ::elem_t
 *
 * @param[in] argc кол-во аргументов (требуется 1)
 * @param[in] argv массив char** содержит {адрес с файлами геометрии}
 * @return int ::e_type_completion
 */
int ReWriteBinToGeo(int argc, char **argv);

/**
 * @brief Функция строит сетку vtk из бинарный файлов
 *
 * @param argc кол-во аргументов (требуется 3)
 * @param argv массив char** содержит {файл vtk сетки, адрес с файлами решения, число сеток}
 * @return int ::e_type_completion
 * @note адрес задаётся с включение главного имени без индекса
 */
int RebuildSolve(int argc, char **argv);

/**
 * @brief Функция строит для всех модулей геометрические настроечные файлы из vtk сетки
 *
 * @param argc кол-во аргументов (требуется 3)
 * @param argv массив char** содержит {адрес для сборки, файл vtk сетки}
 * @return int ::e_type_completion
 */
int MakeAllGeoData(int argc, char **argv);

/**
 * @brief Возвращает средний и минимальный размер ячеек (радиус вписанной окружности)
 *
 * @param argc кол-во аргументов (требуется 1)
 * @param argv массив char** содержит {файл vtk сетки}
 * @return int  ::e_type_completion
 */
int GetAverageSize3D(int argc, char **argv);

/**
 * @brief Создает файл с полной проекцией решения на заданную координатную ось
 *
 * @note у векторных компонент берётся норма
 * @param argc кол-во аргументов (требуется 3)
 * @param argv массив char** содержит {файл vtk сетки, адрес вывода файлов, ось X,Y,Z}
 * @return int ::e_type_completion
 */
int Make1dProjection(int argc, char **argv);

/**
 * @brief Создает серию файлов с полной проекцией решения на заданную координатную ось
 *
 * @note у векторных компонент берётся норма
 * @param argc кол-во аргументов (требуется 4)
 * @param argv массив char** содержит {шаблон названия vtk сетки, адрес вывода файлов, ось X,Y,Z, число проекций(решений)}
 * @return int ::e_type_completion
 */
int MakeSeries1dProjection(int argc, char **argv);
/**
 * @brief Переписывает поля сетки в бинарные файлы
 *
 * @param argc кол-во аргументов (требуется 2)
 * @param argv массив char** содержит {файл vtk сетки, адрес вывода файлов}
 * @return int ::e_type_completion
 */
int WriteBinFromVtk(int argc, char *argv[]);

/**
 * @brief Формирует файл с номерами ячеек, пересекавших луч
 *
 * @param argc кол-во аргументов (требуется 8)
 * @param argv массив char** содержит {файл vtk сетки, файл результат, направление{x,y,z}, основание {x,y,z}}
 * @return int ::e_type_completion
 */
int MakeTrace(int argc, char *argv[]);

/**
 * @brief Создает серию файлов с проекцией (по лучу) решения на заданную координатную ось
 *
 * @param argc кол-во аргументов (требуется 10)
 * @param argv массив char** содержит {шаблон названия vtk сетки, адрес вывода файлов,
 * направление{x,y,z}, основание {x,y,z} число проекций(решений), ось X,Y,Z}
 * @return int ::e_type_completion
 */
int MakeSeries1dRayProjection(int argc, char *argv[]);

/**
 * @brief Формирует файл с гранями всей границы(включая внутреннюю) в формате ::FaceCell
 *
 * @param argc кол-во аргументов (требуется 2)
 * @param argv массив char** содержит {шаблон названия vtk сетки, файл результата}
 * @return int ::e_type_completion
 */
int GetSurface(int argc, char *argv[]);
#endif

} // namespace utils

/**
 * \todo добавить утилиты предсборки, а не функциями
 *
 */

#endif //! UTILS_H
