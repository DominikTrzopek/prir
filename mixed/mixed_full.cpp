#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include <cctype>

#define DEF_MIN_VAL 0
#define DEF_MAX_VAL 1

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

namespace normalization
{
    // Get number of all rows in input file
    int getNumRows(const std::string &filename)
    {
        std::ifstream file(filename);
        int numRows = 0;
        if (file.is_open())
        {
            std::string line;
            while (std::getline(file, line))
            {
                numRows++;
            }
            file.close();
        }
        else
        {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
        return numRows;
    }

    // Load rows from x to y from input file
    std::vector<std::vector<std::string>> loadData(const std::string &filename, int x, int y)
    {
        std::ifstream file(filename);
        std::vector<std::vector<std::string>> data;

        if (file.is_open())
        {
            std::string line;
            int currentRow = 0;

            while (std::getline(file, line))
            {
                if (currentRow < x)
                {
                    // Skip rows until reaching the starting row (x)
                    currentRow++;
                    continue;
                }
                else if (currentRow > y)
                {
                    // Stop reading after reaching the ending row (y)
                    break;
                }

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

                currentRow++;
            }

            file.close();
        }
        else
        {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
        return data;
    }

    // Check if value is numeric
    bool isNumeric(const std::string &str)
    {
        for (char c : str)
        {
            if (!(std::isdigit(c) || c == '.' || c == '-'))
            {
                return false;
            }
        }
        return true;
    }

    // Return min and max value as tuple
    std::vector<double> findMinMaxValues(const std::vector<std::vector<std::string>> &data,  const std::vector<int> &skipCol)
    {
        std::vector<double> result;
        result.reserve(data[0].size() * 2);

        for (size_t col = 0; col < data[0].size(); ++col)
        {
            if(!isNumeric((data[1][col])) || (std::find(skipCol.begin(), skipCol.end(), col) != skipCol.end()))
            {
                result.emplace_back(DEF_MIN_VAL);
                result.emplace_back(DEF_MAX_VAL);
                continue;
            }
            double minVal = std::numeric_limits<double>::max();
            double maxVal = std::numeric_limits<double>::lowest();

            
            #pragma omp parallel for reduction(min:minVal) reduction(max:maxVal)
            for (size_t row = 0; row < data.size(); ++row)
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

            // Handle the case where all values are the same or column is empty
            if (maxVal - minVal < 0.0001)
            {
                minVal = DEF_MIN_VAL;
                maxVal = DEF_MAX_VAL;
            }
            result.emplace_back(minVal);
            result.emplace_back(maxVal);
        }
        return result;
    }

    // Function to normalize data (for numeric values only)
    // Iterate through columns if it contains numeric value convert to double, find min, max, normalize data and save to string
    void normalizeData(std::vector<std::vector<std::string>> &data, std::vector<double> minMaxVec)
    {
        #pragma omp parallel for
        for (size_t col = 0; col < data[0].size(); ++col)
        {
            try
            {
                int it = col * 2;
                double minVal = minMaxVec[it];
                double maxVal = minMaxVec[it+1];

                // Find if column is on the list to not be normalized
                if (minVal == DEF_MIN_VAL && maxVal == DEF_MAX_VAL)
                {
                    continue;
                }

                // Normalize each numeric value in the column
                #pragma omp parallel for
                for (size_t row = 0; row < data.size(); ++row)
                {
                    double value = std::stod(data[row][col]);
                    data[row][col] = std::to_string((value - minVal) / (maxVal - minVal));
                }
            }
            catch (const std::exception &e)
            {
                continue;
            }
        }
    }

    // Function to save data to a CSV file
    void saveData(const std::string &filename, const std::vector<std::vector<std::string>> &data)
    {
        std::ofstream file(filename, std::ios_base::app); // Open the file in append mode

        if (file.is_open())
        {
            for (const auto &row : data)
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


} // namespace normalization

namespace helpers
{
    // Return vector of columns to be skipped during normalization
    std::vector<int> parseColumnsToSkip(int argc, char *argv[])
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
                catch (const std::invalid_argument &e)
                {
                    continue;
                }
            }
        }
        return skipCol;
    }

    // Return min and max value as tuple
    std::pair<int, int> splitRows(int MPI_rank, int MPI_size, int numOfRows)
    {
        int rowsPerProcess = numOfRows / MPI_size;
        int remainingRows = numOfRows % MPI_size;

        int startRow = MPI_rank * rowsPerProcess;
        int endRow = (MPI_rank + 1) * rowsPerProcess - 1;

        // The last process takes care of the remaining rows
        if (MPI_rank == MPI_size - 1)
        {
            endRow += remainingRows;
        }

        return std::make_pair(startRow, endRow);
    }

    // Return vector of global min amd max values
    std::vector<double> reduceReceivedValues(const std::vector<std::vector<double>> &receivedValues)
    {
        std::vector<double> reducedVal;
        for (size_t row = 0; row < receivedValues[0].size(); row += 2)
        {
            double min = std::numeric_limits<double>::max();
            double max = std::numeric_limits<double>::lowest();
            for (size_t col = 0; col < receivedValues.size(); ++col)
            {
                double recMin = receivedValues[col][row];
                double recMax = receivedValues[col][row + 1];

                if(recMin < min)
                {
                    min = recMin;
                }

                if(recMax > max)
                {
                    max = recMax;
                }

            }
            reducedVal.push_back(min);
            reducedVal.push_back(max);
        }
        return reducedVal;
    }

    void clearAndRewriteFirstRow(const std::string &inputFilename, const std::string &outputFilename)
    {
        // Read the first row from the input file
        std::ifstream inputFile(inputFilename);
        if (!inputFile.is_open())
        {
            std::cerr << "Unable to open input file: " << inputFilename << std::endl;
            return;
        }

        std::string firstRow;
        std::getline(inputFile, firstRow);
        inputFile.close();

        // Clear the content of the output file and rewrite the first row
        std::ofstream outputFile(outputFilename, std::ios_base::trunc);
        if (!outputFile.is_open())
        {
            std::cerr << "Unable to open output file: " << outputFilename << " for clearing" << std::endl;
            return;
        }

        outputFile << firstRow << "\n";
        outputFile.close();

        std::cout << "First row from input file rewritten to output file." << std::endl;
    }

} // namespace helpers

// Main function
int main(int argc, char *argv[])
{
    Timer timer;

    int MPI_rank, MPI_size;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

    // Check and parse args
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> -skip <column 1> <column 2> ..." << std::endl;
        return 1;
    }

    std::string inputFilename = argv[1];
    std::string outputFilename = argv[2];
    std::vector<int> skipCol = helpers::parseColumnsToSkip(argc, argv);

    // Load data from CSV file
    int numOfRows = normalization::getNumRows(inputFilename);
    auto [rowsStart, rowsEnd] = helpers::splitRows(MPI_rank, MPI_size, numOfRows);
    if(MPI_rank == 0)
    {
        rowsStart = 1;
    }
    std::vector<std::vector<std::string>> data = normalization::loadData(inputFilename, rowsStart, rowsEnd);

    if (!data.empty())
    {
        // Get min max vector
        std::cout << "P" << MPI_rank << " -> Data loaded from CSV (" << rowsStart << ", " << rowsEnd << ") -> " << timer.printElapsedTime() << " ms " << std::endl;
        std::vector<double> minMaxVec =  normalization::findMinMaxValues(data, skipCol);
        std::vector<double> reducedMinMax(minMaxVec.size());
        
        if(MPI_rank == 0)
        {
            helpers::clearAndRewriteFirstRow(inputFilename, outputFilename);

            std::vector<std::vector<double>> receivedValues;
            receivedValues.push_back(minMaxVec);
            for (int source = 1; source < MPI_size; ++source)
            {
                std::vector<double> buffer(minMaxVec.size());
                MPI_Recv(buffer.data(), minMaxVec.size(), MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                receivedValues.push_back(buffer);
            }

            reducedMinMax = helpers::reduceReceivedValues(receivedValues);
            std::vector<double> reducedMinMaxToSend(reducedMinMax);

            std::cout << "P" << MPI_rank << " -> Min max vector ready -> " << timer.printElapsedTime() << " ms " << std::endl;

            for (int source = 1; source < MPI_size; ++source)
            {
                MPI_Send(reducedMinMaxToSend.data(), reducedMinMaxToSend.size(), MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Send(minMaxVec.data(), minMaxVec.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(reducedMinMax.data(), minMaxVec.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Normalize data
        normalization::normalizeData(data, reducedMinMax);
        std::cout << "P" << MPI_rank << " -> Normalization completed -> " << timer.printElapsedTime() << " ms " << std::endl;

        int ready = 1;
        if(MPI_rank > 0)
        {
            MPI_Recv(&ready, 1, MPI_INT, MPI_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Save to output
        normalization::saveData(outputFilename, data);

        if(MPI_rank < MPI_size - 1)
        {
            MPI_Send(&ready, 1, MPI_INT, MPI_rank + 1, 0, MPI_COMM_WORLD);
        }
        std::cout << "P" << MPI_rank << " -> Data saved to output file -> " << timer.printElapsedTime() << " ms " << std::endl;
    }

    MPI_Finalize();
    return 0;
}