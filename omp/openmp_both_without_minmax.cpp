#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <omp.h>

#define PRECISION 0.0000001

// Time measurement
class Timer
{
public:
    Timer() : start_time(std::chrono::high_resolution_clock::now()) {}

    auto printElapsedTime() const
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        return elapsed_time.count();
    }

private:
    std::chrono::high_resolution_clock::time_point start_time;
};

// Function to load data from a CSV file
// Each column loaded as vector of string values
std::vector<std::vector<std::string>> loadData(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<std::vector<std::string>> data;

    if (file.is_open())
    {
        std::string line;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::vector<std::string> row;
            std::string value;

            while (std::getline(iss, value, ','))
            {
                row.push_back(value);
            }

            if (!row.empty())
            {
                data.push_back(row);
            }
        }

        file.close();
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return data;
}

// Function to normalize data (for numeric values only)
// Iterate through columns if it contains numeric value convert to double, find min, max, normalize data and save to string
void normalizeData(std::vector<std::vector<std::string>>& data, const std::vector<int>& skipCol)
{
    #pragma omp parallel for
    for (size_t col = 0; col < data[0].size(); ++col)
    {
        try
        {
            // Find if column is on the list to not be normalized
            if (std::find(skipCol.begin(), skipCol.end(), col) != skipCol.end())
            {
                throw std::exception();
            }

            // Find min and max values for each column
            double minVal = std::stod(data[1][col]);
            double maxVal = std::stod(data[1][col]);
            for (size_t row = 1; row < data.size(); ++row)
            {
                double value = std::stod(data[row][col]);

                if (value < minVal)
                {
                    minVal = value;
                }

                if (value > maxVal)
                {
                    maxVal = value;
                }
            }

            if (std::abs(maxVal - minVal) < PRECISION) {
                throw std::exception();
            }

            // Normalize each numeric value in the column
            #pragma omp parallel for
            for (size_t row = 1; row < data.size(); ++row)
            {
                double value = std::stod(data[row][col]);
                data[row][col] = std::to_string((value - minVal) / (maxVal - minVal));
            }
        }
        catch (const std::exception& e)
        {
            // Skip non-numeric or user-defined columns
            std::cout << "Skipped column " << col << std::endl;
            continue;
        }
    }
}

// Function to save data to a CSV file
void saveData(const std::string& filename, const std::vector<std::vector<std::string>>& data)
{
    std::ofstream file(filename);

    if (file.is_open())
    {
        for (const auto& row : data)
        {
            for (size_t col = 0; col < row.size(); ++col)
            {
                file << row[col];
                if (col < row.size() - 1)
                {
                    file << ",";
                }
            }
            file << "\n";
        }

        file.close();
    }
    else
    {
        std::cerr << "Unable to create file: " << filename << std::endl;
    }
}

// Return vector of columns to be skipped during normalization
std::vector<int> parseColumnsToSkip(int argc, char* argv[])
{
    std::vector<int> skipCol;
    if (argc > 3 && std::string(argv[3]) == "-skip")
    {
        for (int i = 4; i < argc; ++i)
        {
            try
            {
                skipCol.push_back(std::stoi(argv[i]));
            }
            catch (const std::invalid_argument& e)
            {
                continue;
            }
        }
    }
    return skipCol;
}

// Main function
int main(int argc, char* argv[])
{
    Timer timer;

    // Check and parse args
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> -skip <column 1> <column 2> ..." << std::endl;
        return 1;
    }
    std::string inputFilename = argv[1];
    std::string outputFilename = argv[2];
    std::vector<int> skipCol = parseColumnsToSkip(argc, argv);

    std::cout << "Program starting -> " << timer.printElapsedTime() << " ms" << std::endl;

    // Load data from CSV file
    std::vector<std::vector<std::string>> data = loadData(inputFilename);

    if (!data.empty())
    {
        std::cout << "Data loaded from CSV -> " << timer.printElapsedTime() << " ms" << std::endl;

        // Normalize data
        normalizeData(data, skipCol);
        std::cout << "Normalization completed -> " << timer.printElapsedTime() << " ms" << std::endl;

        // Save data to CSV file
        saveData(outputFilename, data);
        std::cout << "Output saved to " << outputFilename << " -> " << timer.printElapsedTime() << " ms" << std::endl;
    }

    return 0;
}
