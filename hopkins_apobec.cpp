#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <time.h>
#include <cmath>
#include <thread>
using namespace std;

void wholeProgram( string chr,  int idIndex, vector<string>  idList, vector< vector<string> >  data, std::ofstream & output);
vector<int> generateRandomIntegers(int n, int R);
vector< vector<string> > read_tsv(string directory);
vector<int> deleteClusters(vector<int> & raw);
bool compare(const vector<string>& a, const vector<string>& b);
vector<string> getRowOfGene(int position, const vector< vector<string> > & data);
struct coordination;
double getDistanceEuclidean3D(coordination coordination1, coordination coordination2);

vector<int> generateRandomIntegers(int n, int R) 
{
	// Create a random number generator
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(0, R);

	// Generate n random integers between 0 and R
	vector<int> result(n);
	for (int i = 0; i < n; i++) {
		result[i] = dis(gen);
	}

	// Return the vector of random integers
	return result;
}

vector< vector<string> > read_tsv(string directory)
{
    string line;
    ifstream myfile(directory);
    vector< vector<string> > result;
    if (myfile.is_open())
    {
        while (getline(myfile,line))
        {
            vector<string> row;
            int start = 0;
            for(int i = 0; i<line.size(); i++)
            {
                if(line[i] == '\t')
                {
                    row.push_back(line.substr(start, i-start));
                    start = i+1;
                }
                else if((i + 1) == line.size())
                {
                    row.push_back(line.substr(start, i+1-start));
                }
            }
            result.push_back(row);
        }
        myfile.close();
    }
    else cout << "Unable to open file"; 

    return result;
}

vector<int> deleteClusters(vector<int> & raw)
{
    int lastMutation = 0;
    int currentClusterSize = 0;
    vector<int> result;

    int prevEnd = 0;
    int currentStart;
 
    for(int i = 0; i < raw.size(); i++)
    {
        if(((raw[i] - lastMutation) >= 10000) && (lastMutation != 0))
        {
            if(currentClusterSize >= 2)
            {
                currentStart = i - currentClusterSize;
                for(int a = (prevEnd+1); a < currentStart; a++)
                {
                    result.push_back(raw[a]);
                }
                prevEnd = i - 1;
            }
            currentClusterSize = 1;
            lastMutation = raw[i];
        }
        else
        {
            currentClusterSize++;
            lastMutation = raw[i];
        }
    }
    return result;
}

bool compare(const vector<string>& a, const vector<string>& b)
{
    return stoi(a[3]) < stoi(b[3]); // Compare the 4th value of each vector
}

vector<string> getRowOfGene(int position, const vector< vector<string> > & data) 
{
    int low = 0;
    int high = data.size() - 1;

    while (low <= high) 
    {
        int mid = low + (high - low) / 2;
        int start = stoi(data[mid][0]);
        int end = stoi(data[mid][1]);
        if (position >= start && position <= end)return data[mid];
        else if (position < start) high = mid - 1;
        else low = mid + 1;
    }

    return vector<string>();
}

struct coordination
{
    double x;
    double y;
    double z;
    bool empty;
    coordination(int position, string chr, vector< vector<string> > hsa)
    {
        if(getRowOfGene(position,  hsa).size() != 0)
        {
            x = stod(getRowOfGene(position,  hsa)[2]);
            y = stod(getRowOfGene(position,  hsa)[3]);
            z = stod(getRowOfGene(position,  hsa)[4]);
            empty = false;
        }
        else
        {
            x = 0.0;
            y = 0.0;
            z = 0.0;
            empty = true;
        }
    };
};


double getDistanceEuclidean3D(coordination coordination1, coordination coordination2)
{
    double x = (coordination1.x - coordination2.x)*(coordination1.x - coordination2.x);
    double y = (coordination1.y - coordination2.y)*(coordination1.y - coordination2.y);
    double z = (coordination1.z - coordination2.z)*(coordination1.z - coordination2.z);
    return sqrt(x+y+z);
}

void wholeProgram(string chr, int idIndex, vector<string> idList, vector< vector<string> > data, std::ofstream & output)
{

    string id = idList[idIndex];
	vector< vector<string> > target_ID_Chr;
	for(int i = 1; i<data.size(); i++)
	{
        if(data[i][1] == id)
        {
            if(data[i][2] == chr)
                target_ID_Chr.push_back(data[i]);
        }
	}
    if (target_ID_Chr.empty())
    {
        output << idList[idIndex] <<"\t" << chr <<"\t"<< "-" <<"\t"<<  "0"  <<"\t"<<  "0" <<"\t"<< "0" << endl;
        return;
    }
	vector<int> indexOfMutations;
	for(int i = 0; i<target_ID_Chr.size(); i++)
	{
		if(target_ID_Chr[i][8] == "1")
			indexOfMutations.push_back(i);
	}
	
    if (indexOfMutations.empty())
    {
        output << idList[idIndex] <<"\t" << chr <<"\t"<< "-" <<"\t"<<  0  <<"\t"<<  "0" <<"\t"<< "0" << endl;
        return;
    }
    vector<int> allMutations;
    for(int i = 0; i < indexOfMutations.size()-1; i++)
    {
        allMutations.push_back(stoi(target_ID_Chr[indexOfMutations[i]][3]));
    }

    
    allMutations.push_back(stoi(target_ID_Chr[indexOfMutations[(int)indexOfMutations.size()-1]][3]));


    vector<int> scattereds = deleteClusters(allMutations);
    // cout << "Count of scattered mutations == " << scattereds.size() <<"."<<endl;

    int sampleSize =  (((int)scattereds.size())/20);
    if(sampleSize < 100) sampleSize = 100;
    if(((int)scattereds.size()) < sampleSize) sampleSize = ((int)scattereds.size());
    // cout << "Sample size for hopkins statistics is   " << sampleSize << "." << endl;

    if(scattereds.size() < sampleSize) sampleSize = (int) scattereds.size();
    double sum = 0;
    int count = 1;


    for(int k = 0; k<count; k++)
    {
        vector<int> sampleIndexs = generateRandomIntegers(sampleSize, (int) scattereds.size());
        sort(sampleIndexs.begin(), sampleIndexs.end());
        vector<int> sample;
        for(int i = 0; i < sampleIndexs.size(); i++)
        {
            sample.push_back(scattereds[sampleIndexs[i]]);
        }

        string directory2 = "/Users/melih/Desktop/Hi-C_from_A549/HSA/";
        vector< vector<string> > hsa = read_tsv(directory2 + "chr" + chr + "_HSA_out.txt");
        
        int MAX_COORD = stoi(hsa[hsa.size()-1][1]);
        //cout << "MAX_COORD: " << MAX_COORD << endl;
        vector<int> uniform = generateRandomIntegers(sampleSize, MAX_COORD);
        
        double u = 0;
        for (double i = 0; i < uniform.size(); i++)
        {
            coordination uniformC = coordination(uniform[i], chr, hsa);
            while(uniformC.empty)
            {
                random_device rd;
	            mt19937 gen(rd());
	            uniform_int_distribution<> dis(0, MAX_COORD);
                uniformC = coordination(dis(gen), chr, hsa);
            }
            double min = 10000000000.0;
            for(int c=0; c < scattereds.size(); c++)
            {
                coordination mutC = coordination(scattereds[c], chr, hsa);
                if (mutC.empty) continue;
                double dist = getDistanceEuclidean3D(uniformC, mutC);
                if(dist < min) min = dist;
            }
            /* cout << (i/(2*sample.size()))*100 << "% is completed." << endl;
            clock_t end_time = clock();
            double duration_sec = double(end_time - start_time) / CLOCKS_PER_SEC; // calculate duration in seconds
            cout << "Program duration: " << duration_sec << " seconds." << endl; 
            cout << "miminimum distance  '" << i << "' (random):    "<< min << endl;  */
 
            u = u + (min*min*min);
        }
            
        double w = 0;
        for (double i = 0; i < sample.size(); i++)
        {
            coordination sampleC = coordination(sample[i], chr, hsa);
            if(sampleC.empty) continue;
            double min = 1000000000000.0;
            for(int c=0; c < scattereds.size(); c++)
            {
                if(scattereds[c] == sample[i]) continue;
                coordination mutC = coordination(scattereds[c], chr, hsa);
                if (mutC.empty) continue;
                double dist = getDistanceEuclidean3D(sampleC, mutC);
                if(dist < min) min = dist;
            }
            /* cout << (0.5 + (i/(2*sample.size())))*100 << "% is completed." << endl;
            clock_t end_time = clock();
            double duration_sec = double(end_time - start_time) / CLOCKS_PER_SEC; // calculate duration in seconds
            cout << "Program duration: " << duration_sec << " seconds." << endl; 
            cout << "miminimum distance  '" << i << "' (real):    "<< min << endl; */

            w = w + (min*min*min);
        }
        double hopkins = u/(u+w);
        //cout << "Iteration number: " << (k+1) << ". Hopkins Statistics is:     " << hopkins << "." << endl << endl << endl;
        sum += hopkins;

        output << idList[idIndex] <<"\t" << chr <<"\t"<< hopkins <<"\t"<<  allMutations.size()  <<"\t"<<  scattereds.size() <<"\t"<< sampleSize << endl;
    

    }

    
    
}


int main()
{
    clock_t start_time = clock();
    string file;
    string directory = "/Users/melih/desktop/APOBECflag/";
    //cout << "Input file name: ";
    //cin >> file;
    file = "snv_mnv_LUAD-US_short.tsv";
    //cout << endl;

    vector< vector<string> > data = read_tsv(directory + file);
    sort(data.begin()+1, data.end(), compare);
//******************************************

    vector<string> idList;
    for(int i = 1; i<data.size(); i++)
    {
        bool found = false;
        int a = 0;
        while (a < idList.size())
        {
            if(data[i][1] == idList[a])
            {
                found = true;
            }
            a++;
        }
        if(!found) 
        {
            idList.push_back(data[i][1]);
            //cout << data[i][1] << endl;
        }
    }

    std::ofstream finalOutput("finalOutput.txt", std::ios::app);
    finalOutput << "ID" << "\t" << "CHR" << "\t" << "HOPKINS" << "\t" << "MUTATION_COUNT"<< "\t" << "SCATTERED_MUTATION_COUNT"<< "\t" << "SAMPLE_SIZE" << endl;
    
    std::vector<std::thread> vecThreads;

    // create and launch the threads
    for (int i = 0; i < idList.size(); ++i) {
        for (int a = 1; a < 23; a++)
        vecThreads.emplace_back(wholeProgram, to_string(a), i, idList, data, ref(finalOutput));
        string x = "X";
        string y = "Y";
        vecThreads.emplace_back(wholeProgram, x, i, idList, data, ref(finalOutput));
        vecThreads.emplace_back(wholeProgram, y, i, idList, data, ref(finalOutput));

        // wait for the threads to finish
        for (auto& thread : vecThreads) {
            thread.join();
        }
        vecThreads.clear();
    }

    

    finalOutput.close();
    
    clock_t end_time = clock();
    double duration_sec = double(end_time - start_time) / CLOCKS_PER_SEC; // calculate duration in seconds
    cout << "Program duration: " << duration_sec << " seconds." << endl;
    
    return 0;
}
