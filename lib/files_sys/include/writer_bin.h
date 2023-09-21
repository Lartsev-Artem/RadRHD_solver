#if 0
#define WRITE_FILE_VECTOR(name_file, data, value)   \
  {                                                 \
    FILE *f;                                        \
    int n = data.size();                            \
    f = fopen(name_file, "wb");                     \
    if (!f)                                         \
      RETURN_ERRS("file %s not open\n", name_file); \
    fwrite(&n, sizeof(int), 1, f);                  \
    for (auto &el : data) {                         \
      fwrite(&el.value, sizeof(el.value), 1, f);    \
    }                                               \
    fclose(f);                                      \
  }

#define WRITE_FILE(name_file, data, n)              \
  {                                                 \
    FILE *f;                                        \
    f = fopen(name_file, "wb");                     \
    if (!f)                                         \
      RETURN_ERRS("file %s not open\n", name_file); \
    fwrite(&n, sizeof(int), 1, f);                  \
    fwrite(data, sizeof(data[0]), n, f);            \
    fclose(f);                                      \
  }

#define WRITE_FILE_PHYS(name_file, data, value, type, param) \
  {                                                          \
    FILE *f;                                                 \
    int n = data.size();                                     \
    f = fopen(name_file, "wb");                              \
    if (!f)                                                  \
      RETURN_ERRS("file %s not open\n", name_file);          \
    fwrite(&n, sizeof(int), 1, f);                           \
    for (auto &el : data) {                                  \
      type x = el.value * param;                             \
      fwrite(&x, sizeof(x), 1, f);                           \
    }                                                        \
    fclose(f);                                               \
  }
#endif