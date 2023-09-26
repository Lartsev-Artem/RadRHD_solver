#include "global_def.h"
#include "utils.h"
#include <iterator>

/// возможно стоит упростить до чтения из файла в файл
//   int buf;
//   int n = 0;
//   std::list<int> list;
//   while (ifile >> buf) {
//     list.add(buf);
//     n++;
//   }
//   for (auto el : list) {
//     ofile << el << '\n';
//   }

int FUNC_NAME(DataToSimpleTxt)(int argc, char **argv) {
  if (argc != 3) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\data.txt\n");
    printf("path\\out_simple_format.txt\n");
    return 1;
  }

  std::string file_data = argv[1];
  std::string file_simple_format = argv[2];

  // std::string file_data = "D:\\Desktop\\FilesCourse\\Test\\auto_metis\\cone_m_epart_2.txt";
  // std::string file_simple_format = "D:\\Desktop\\FilesCourse\\Test\\auto_metis\\cone_m_epart_2_simple.txt";

  std::ifstream ifile;
  OPEN_FSTREAM(ifile, file_data.c_str());

  std::string file_content; // ((std::istream_iterator<char>(ifile.rdbuf())), std::istream_iterator<char>());
  file_content.assign((std::istreambuf_iterator<char>(ifile.rdbuf())), std::istreambuf_iterator<char>());
  ifile.close();

  OPEN_FSTREAM(ifile, file_data.c_str());
  int n = std::count(std::istreambuf_iterator<char>(ifile), std::istreambuf_iterator<char>(), '\n');
  ifile.close();

  file_content = std::to_string(n) + "\n" + file_content;

  std::ofstream ofile;
  OPEN_FSTREAM(ofile, file_simple_format.c_str());
  ofile << file_content;
  ofile.close();

  return e_completion_success;
}
