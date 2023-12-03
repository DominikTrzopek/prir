#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include <cctype>

#define INVALID_COLUMN -1

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
    // Function to load data from a CSV file
    // Each column loaded as vector of string values
    std::vector<std::vector<std::string>> loadColumns(const std::string &filename, const std::vector<int> &columnIndexes)
    {
        std::ifstream file(filename);
        std::vector<std::vector<std::string>> columnsData;

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

                // Check if the row is not empty
                if (!row.empty())
                {
                    // Extract columns at specified indexes
                    std::vector<std::string> selectedColumns;
                    for (int index : columnIndexes)
                    {
                        if (std::abs(index) < row.size())
                        {
                            selectedColumns.push_back(row[std::abs(index)]);
                        }
                        else
                        {
                            // Handle the case where the index is out of bounds for the current row
                            selectedColumns.push_back("");
                        }
                    }

                    columnsData.push_back(selectedColumns);
                }
            }

            file.close();
        }
        else
        {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }

        return columnsData;
    }

    // Return min and max value as tuple
    std::pair<double, double> findMinMaxValues(const std::vector<std::vector<std::string>> &data, int col)
    {
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

        if (maxVal - minVal < 0.0001)
        {
            minVal = 0;
            maxVal = 1;
        }

        return std::make_pair(minVal, maxVal);
    }

    // Function to normalize data (for numeric values only)
    // Iterate through columns if it contains numeric value convert to double, find min, max, normalize data, save as double
    std::vector<std::vector<std::string>> normalizeData(const std::vector<std::vector<std::string>> &data)
    {
        std::vector<std::vector<std::string>> normalized(data[0].size(), std::vector<std::string>(data.size(), ""));

        for (size_t col = 0; col < data[0].size(); ++col)
        {
            try
            {
                // Find min and max values for each column
                auto [minVal, maxVal] = findMinMaxValues(data, col);

                // Normalize each numeric value in the column
                for (size_t row = 1; row < data.size(); ++row)
                {
                    double value = std::stod(data[row][col]);
                    normalized[col][row] = std::to_string((value - minVal) / (maxVal - minVal));
                }
            }
            catch (const std::exception &e)
            {
                for (size_t row = 1; row < data.size(); ++row)
                {
                    normalized[col][row] = data[row][col];
                }
                continue;
            }
        }

        return normalized;
    }

    // Function to normalize data (for numeric values only)
    // Iterate through columns if it contains numeric value convert to double, find min, max, normalize data, save as double
    std::vector<std::vector<double>> normalizeDataDouble(const std::vector<std::vector<std::string>> &data)
    {
        std::vector<std::vector<double>> normalized(data[0].size(), std::vector<double>(data.size(), 0.0));

        for (size_t col = 0; col < data[0].size(); ++col)
        {
            // Find min and max values for each column
            auto [minVal, maxVal] = findMinMaxValues(data, col);

            // Normalize each numeric value in the column
            for (size_t row = 1; row < data.size(); ++row)
            {
                double value = std::stod(data[row][col]);
                normalized[col][row] = (value - minVal) / (maxVal - minVal);
            }
        }

        return normalized;
    }

    // Save data to csv, first save column names from input file, then write data
    void saveData(const std::string &filename, const std::string &inputFilename, const std::vector<std::vector<std::string>> &data)
    {
        std::ofstream file(filename);

        if (file.is_open())
        {
            // Get the number of rows and columns
            size_t numRows = data.size();
            size_t numCols = (numRows > 0) ? data[0].size() : 0;

            // Check if there is at least one row in the input data
            if (numRows > 0)
            {
                // Open the input file
                std::ifstream inputFile(inputFilename);

                if (inputFile.is_open())
                {
                    // Read the first row from the input file
                    std::string firstRow;
                    std::getline(inputFile, firstRow);

                    // Write the first row from the input file to the output file
                    file << firstRow << "\n";

                    // Iterate over columns first, then rows starting from the second row
                    for (size_t col = 1; col < numCols; ++col)
                    {
                        for (size_t row = 0; row < numRows; ++row)
                        {
                            file << data[row][col];
                            if (row < numRows - 1)
                            {
                                file << ",";
                            }
                        }
                        file << "\n";
                    }

                    inputFile.close();
                }
                else
                {
                    std::cerr << "Unable to open input file: " << inputFilename << std::endl;
                }
            }
            else
            {
                std::cerr << "Input data is empty." << std::endl;
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

    // Return second row od data
    std::vector<std::string> getSecondRow(const std::string &csvFile)
    {
        try
        {
            std::ifstream file(csvFile);
            if (!file.is_open())
            {
                std::cerr << "Error: File '" << csvFile << "' not found." << std::endl;
                return {};
            }

            std::string line;
            std::vector<std::string> secondRow;

            // Skip the first row
            if (std::getline(file, line))
            {
                // Read the second row
                if (std::getline(file, line))
                {
                    std::istringstream stream(line);
                    std::string cell;

                    while (std::getline(stream, cell, ','))
                    {
                        secondRow.push_back(cell);
                    }
                }
                else
                {
                    // No second row found
                    std::cerr << "Error: No second row in file '" << csvFile << "'." << std::endl;
                }
            }
            else
            {
                // No first row found
                std::cerr << "Error: No rows in file '" << csvFile << "'." << std::endl;
            }

            file.close();
            return secondRow;
        }
        catch (const std::exception &e)
        {
            std::cerr << "An error occurred: " << e.what() << std::endl;
            return {};
        }
    }

    // Check if string contains non numeric values
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

    // Return vector of columns that can be normalized
    std::vector<int> getValidColumns(int argc, char *argv[], const std::string &csvFile)
    {
        std::vector<int> skipCol = parseColumnsToSkip(argc, argv);
        std::vector<std::string> secondRow = getSecondRow(csvFile);
        int col = 0;

        std::vector<int> validColumns{};

        for (std::string input : secondRow)
        {
            if ((std::find(skipCol.begin(), skipCol.end(), col) != skipCol.end()) || !isNumeric(input))
            {
                // Push back -1 to indicate column should be skipped during normalization
                validColumns.push_back(INVALID_COLUMN);
                col++;
                continue;
            }

            validColumns.push_back(col);
            col++;
        }

        return validColumns;
    }

    // Return vector of columns to normalize for each process
    std::vector<std::vector<int>> splitColumnsForProcesses(std::vector<int> &columns, const int numOfProcesses)
    {
        std::vector<std::vector<int>> columnSplit(numOfProcesses);
        // Process 0 receives all columns that should be skipped in order to copy them into the output file
        for (auto it = columns.begin(); it != columns.end(); ++it)
        {
            if (*it == INVALID_COLUMN)
            {
                columnSplit[0].push_back(std::distance(columns.begin(), it));
            }
        }

        // Delete from vector of columns all invalid elements
        columns.erase(std::remove(columns.begin(), columns.end(), INVALID_COLUMN), columns.end());

        int numOfColumns = columns.size();
        int colPerProcess = numOfColumns / (numOfProcesses - 1);

        int iter = 0;
        int procIter = 0;
        int process = 1;
        // Split columns for each process
        while (iter < numOfColumns)
        {
            if (procIter >= colPerProcess)
            {
                process++;
                if (process > numOfProcesses - 1)
                    break;
                procIter = 0;
            }
            columnSplit[process].push_back(columns.at(iter));
            iter++;
            procIter++;
        }

        // Process 0 gets the remaining columns
        columnSplit[0].insert(columnSplit[0].end(), columns.begin() + iter, columns.end());

        return columnSplit;
    }

    // Print indexes of columns assign to each process
    void printDataSplit(const std::vector<std::vector<int>> &split)
    {
        int pro = 0;
        for (const auto &innerVector : split)
        {
            std::cout << "P" << pro << " -> ";
            for (int value : innerVector)
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
            pro++;
        }
    }

    // Verify if every process will receive data to normalize
    bool isDataSplitCorrect(const std::vector<std::vector<int>> &split)
    {
        for (std::vector<int> vec : split)
        {
            if (vec.size() == 0)
                return false;
        }
        return true;
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

    int numOfAllColumns;
    int numOfColumnsSplit;
    int splitSizeForProcess;
    std::vector<int> columns{};

    if (MPI_rank == 0)
    {
        std::cout << "Program starting -> " << timer.printElapsedTime() << " ms" << std::endl;

        // Prepare data split (indexes of columns to be normalized by each process)
        columns = helpers::getValidColumns(argc, argv, inputFilename);
        numOfAllColumns = columns.size();
        std::vector<std::vector<int>> split = helpers::splitColumnsForProcesses(columns, MPI_size);
        if (!helpers::isDataSplitCorrect(split))
        {
            std::cerr << "WRN: To many processes! Some of them won't receive data" << std::endl;
        }
        std::cout << "Data split ready -> " << timer.printElapsedTime() << " ms" << std::endl;
        helpers::printDataSplit(split);
        splitSizeForProcess = split.at(1).size();

        // Send indexes of columns that each process should normalize
        for (int dest_rank = 1; dest_rank < MPI_size; ++dest_rank)
        {
            // Send number of columns to normalize
            MPI_Send(&splitSizeForProcess, 1, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);

            // Send indexes of columns
            MPI_Send(split.at(dest_rank).data(), splitSizeForProcess, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);
        }

        numOfColumnsSplit = split.at(0).size();
        columns.resize(numOfColumnsSplit);
        columns = split.at(0);
    }
    else
    {
        // Receive number of columns to normalize
        MPI_Recv(&numOfColumnsSplit, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Receive indexes of columns to normalize
        columns.resize(numOfColumnsSplit);
        MPI_Recv(columns.data(), numOfColumnsSplit, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // std::cout << "P" << MPI_rank << " -> Received indexes from process 0 -> " << timer.printElapsedTime() << " ms " << std::endl;
    }

    // Load data from CSV file
    std::vector<std::vector<std::string>> data = normalization::loadColumns(inputFilename, columns);

    if (!data.empty())
    {
        std::cout << "P" << MPI_rank << " -> Data loaded from CSV -> " << timer.printElapsedTime() << " ms " << std::endl;

        if (MPI_rank > 0)
        {

            std::vector<std::vector<double>> normalized = normalization::normalizeDataDouble(data);
            std::cout << "P" << MPI_rank << " -> Normalization completed -> " << timer.printElapsedTime() << " ms" << std::endl;

            // Send normalized data to process 0
            for (int i = 0; i < numOfColumnsSplit; i++)
            {
                // Send column index
                MPI_Send(&columns.at(i), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

                // Serialize the vector of strings
                std::vector<double> toSend = normalized.at(i);

                int vectorSize = toSend.size();
                // Send the size of the serialized vector
                MPI_Send(&vectorSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

                // Send the serialized vector
                MPI_Send(toSend.data(), vectorSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            // std::cout << "P" << MPI_rank << " -> Sent normalized data " << timer.printElapsedTime() << " ms" << std::endl;
        }
        else if (MPI_rank == 0)
        {

            std::vector<std::vector<std::string>> normalized = normalization::normalizeData(data);
            std::cout << "P" << MPI_rank << " -> Normalization completed -> " << timer.printElapsedTime() << " ms" << std::endl;

            std::vector<std::vector<std::string>> normalizedDataRec(numOfAllColumns);
            int colNumber;

            for (int source_rank = 1; source_rank < MPI_size; source_rank++)
            {
                for (int i = 0; i < splitSizeForProcess; i++)
                {
                    // Receive column index
                    MPI_Recv(&colNumber, 1, MPI_INT, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    // Receive vector size
                    int vectorSize;
                    MPI_Recv(&vectorSize, 1, MPI_INT, source_rank, 0, MPI_COMM_WORLD, &status);

                    // Allocate memory to receive the vector
                    std::vector<double> recVector(vectorSize);

                    // Receive the vector
                    MPI_Recv(recVector.data(), vectorSize, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, &status);

                    std::vector<std::string> recStringVector(vectorSize);
                    std::transform(recVector.begin(), recVector.end(), recStringVector.begin(), [](double value)
                                   { return std::to_string(value); });

                    // Merge vector with rest of data
                    normalizedDataRec[colNumber].insert(normalizedDataRec[colNumber].end(), recStringVector.begin(), recStringVector.end());
                }
            }

            std::cout << "All data received -> " << timer.printElapsedTime() << " ms" << std::endl;

            // Merge own data (from process 0)
            for (size_t i = 0; i < columns.size(); ++i)
            {
                int colNumber = columns.at(i);
                std::vector<std::string> stringsToInsert = normalized.at(i);

                // Ensure that normalizedData has enough space
                if (colNumber >= normalizedDataRec.size())
                {
                    normalizedDataRec.resize(colNumber + 1);
                }

                normalizedDataRec[colNumber].insert(normalizedDataRec[colNumber].end(), stringsToInsert.begin(), stringsToInsert.end());
            }

            // Save to output file
            normalization::saveData(outputFilename, inputFilename, normalizedDataRec);
            std::cout << "Output saved to " << outputFilename << " -> " << timer.printElapsedTime() << " ms" << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
