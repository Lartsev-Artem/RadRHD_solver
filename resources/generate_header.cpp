#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include <string>
#include <vector>
#include <map>

#include <sstream>

#include <math.h>

#define FIND_WORD(st, word){ \
if((int)st.str().find(word)>=0){\
std::string _val; \
 do \
  {\
    st >> _val;\
  } while (_val != word && _val != ";");\
}\
}

#include <filesystem>
namespace fs = std::filesystem;

bool is_number(const std:: string& str)
{
    return str.find_first_not_of(".0123456789") == std::string::npos;
}

double conv_value(std::string& word, const std::map<std::string, double>& name_val)
{
    double value=0;
    if (is_number(word))
    {
        value = std::stod(word);
    }
    else
    {
        auto pos = name_val.find(word);
        if (pos != name_val.end())
        {
            value = pos->second;
        }
        else
        {
            abort();
        }
    }
    return value;
}

double get_value(std::string& str, const  std::map<std::string, double>& name_val)
{

    const std::string separators{ " /*+-" }; // разделители слов    
    size_t start = str.find_first_of(" ") + 1;
    double value = 1;
    bool is_first = true;
    char c = ' ';
    while (start != std::string::npos) // если нашли слово
    {
        size_t end = str.find_first_of(separators, start + 1); // находим, где кончается слово        

        if (end == std::string::npos)
        {
            end = str.length();
        }

        std::string word = str.substr(start, end-1);

        size_t pos2 = word.find_first_of("/*+-", 0);        
        if (pos2 != 0)
        {
            word = word.substr(0, pos2-1);
        }        

        if (is_first)
        {            
            value = conv_value(word, name_val);          
            is_first = false;
        }
        else
        {            
            switch (c)
            {
            case '*':
                value *= conv_value(word, name_val);
                break;
            case '/':
                value /= conv_value(word, name_val);
                break;
            case '+':
                value += conv_value(word, name_val);
                break;
            case '-':
                value -= conv_value(word, name_val);
                break;
            default:
                abort();
            }
        }

        int pos = str.find_first_of("/*+-", end);
        if (pos >=0)
        {
            c = str[pos];
        }

      
        start = str.find_first_not_of(separators, end + 1); // находим начальный индекс следующего слова и переустанавливаем start
    }
  
    return value;
}

int main(int argc, char** argv)
{

    if (argc != 3)
    {
        printf("Error input format:\n");
        printf("input: input_file, out_file\n");
        return 1;
    }

    const std::string input = argv[1];
    const std::string output = argv[2];

    //const std::string input = "plank.txt";
    //const std::string output = "test.h";

    std::ifstream ifile;
    ifile.open(input);

    std::vector< std::string> types;
    std::vector< std::string> names;
    std::vector<double> values;
    std::map<std::string, double> name_val;
    
    while (!ifile.eof())
    {
        std::string str;
        std::getline(ifile, str);

        std::stringstream ss(str.substr(0, str.find(";") + 1));
        std::string type, name;
        FIND_WORD(ss, "constexpr");
        ss >> type;
        ss >> name;

        bool flag=false;
        double val=0;
        {
            std::string buf = str.substr(str.find("=") + 1, str.find(";") -str.find("=")-1);
            
            for (auto &el: name_val)
            {
                if (buf.find(el.first) != std::string::npos)
                {
                    val= get_value(buf, name_val);
                    flag = true;
                    break;
                }
            }
            
        }

        if (!flag)
        {
            double val2 = 0;
            FIND_WORD(ss, "=");
            if (!(ss >> val))
            {
                continue;
            }

            FIND_WORD(ss, "*");
            if (ss >> val2)
            {
                val *= val2;
            }
        }

        types.push_back(type);
        names.push_back(name);
        values.push_back(val);   
        name_val.insert(std::make_pair(name, val));        
    }
    ifile.close();
    //========================================

    std::ofstream ofile;
    ofile.open(output);

    std::string def_hed = fs::path(output).stem().string(); //output.substr(0, output.find("."));
    std::transform(def_hed.begin(), def_hed.end(), def_hed.begin(), ::toupper);

    ofile << "#ifndef " << def_hed << "_H\n";
    ofile << "#define " << def_hed << "_H\n";

    for (size_t i = 0; i < types.size(); i++)
    {
        ofile << std::setprecision(16) << "constexpr " << types[i] << " log_" << names[i] << " = " << log(values[i]) << ";\n";
        std::cout << std::setprecision(16) << "constexpr " << types[i] << " log_" << names[i] << " = " << values[i] << "->" << log(values[i]) << "; \n";
    }
    ofile << "#endif //! " << def_hed << "_H\n";

    ofile.close();
    return 0;
}


