/**
 * @file file_sys.h
 * @brief Файл подключает все заголовки файлового модуля
 * @version 0.1
 * @date 2023-09-25
 *
 */
#ifndef FILE_SYS_H
#define FILE_SYS_H

#include "reader_bin.h"
#include "reader_json.h"
#include "reader_txt.h"
#include "reader_vtk.h"

#include "writer_bin.h"
#include "writer_json.h"
#include "writer_txt.h"
#include "writer_vtk.h"

/*! \addtogroup file_sys Файловый модуль
    \brief Модуль содержит функции работы с  файлами в различных форматах
    \note поддерживаются txt, bin, json, vtk форматы
    @{
*/

/**
 * @brief Пространство имён файлового модуля
 *
 */
namespace files_sys {}

#endif //! FILE_SYS_H