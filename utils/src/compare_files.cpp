#include "global_def.h"

#include "utils.h"

#include <vector>

#include "geo_types.h"
#include "reader_bin.h"

#include <limits>

typedef double com_type;

int FUNC_NAME(CompareFiles)(int argc, char *argv[]) {
  std::string name_file_1 = "";
  std::string name_file_2 = "";

  if (argc <= 2) {
    RETURN_ERR("Error input data!\n Input: path\\name_file_1, path\\name_file_2\n");
  } else {
    name_file_1 = argv[1];
    name_file_2 = argv[2];
  }

  std::vector<com_type> arr1;
  std::vector<com_type> arr2;
  files_sys::bin::ReadSimple(name_file_1, arr1);
  files_sys::bin::ReadSimple(name_file_2, arr2);

  if (arr1.size() != arr2.size()) {
    RETURN_ERR("Error arrays is different\n");
  }

  com_type max_err = 0;
  int max_idx = 0;

  com_type min_value = std::numeric_limits<com_type>::max();
  com_type max_value = -std::numeric_limits<com_type>::max();
  for (size_t i = 0; i < arr1.size(); i++) {

    if (std::isnan(arr1[i]) || std::isnan(arr2[i]) ||
        std::isinf(arr1[i]) || std::isinf(arr2[i])) {
      std::cout << "Error value[" << i
                << "]: " << arr1[i]
                << " " << arr2[i] << "\n";
      return e_completion_fail;
    }

    min_value = std::min(min_value, std::min(arr1[i], arr2[i]));
    max_value = std::max(max_value, std::max(arr1[i], arr2[i]));

    com_type tmp = fabs(arr1[i] - arr2[i]);
    if (tmp > max_err) {
      max_err = tmp;
      max_idx = i;
    }
  }

  std::cout << "Max err[" << max_idx << "]= " << max_err << "\n";
  std::cout << "Eps err= " << max_err / min_value << "\n";
  std::cout << "values: " << arr1[max_idx] << " " << arr2[max_idx] << "\n";
  std::cout << "MinMax: " << min_value << " " << max_value << "\n";

  return e_completion_success;
}
