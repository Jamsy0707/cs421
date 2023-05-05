/*
    James Sieben 200455325
    CS421 A2 8puzzle solved using the A* algorithm.
    
    To use, simply run an enter an initial puzzle state using no spaces and a 0 for the empty tile.
    ex. 123456780 is the goal state

    Based on code by user Arty on stackoverflow.
*/

#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using u8 = uint8_t;                 //String literal

size_t constexpr n = 3, m = 3;      //Defines puzzle size.                               

int Solve(string const &, string const &);

//The main driver function.
int main() 
{
    string initialArr;              //The puzzle to be solved.
    string finalArr = "123456780";  //The goal state.
    int numStates;

    cout << "This is an implementation of the 8-puzzle problem with the \"Tiles out of place\" heuristic" << endl
         << "Please enter an initial state (Enter 0 for the blank spot, no spaces!):" << endl;

    cin >> initialArr;

    cout << "The sequence is:\n";

    numStates = Solve(initialArr, finalArr);    //Call the solve function using the defined initial and goal states.

    cout << "The number of states is " << numStates << ".\n";

    return 0;
}

/*
  Parameters: N/A
  Returns: N/A
  Description: This class defines functions related to the puzzle boards.
*/
class Puzzle {
    public:
        Puzzle() : n_(n), m_(m), puzzle_(n_ * m_) {} //Define the puzzle size based on the variables declared earlier

        Puzzle(string const & s) : n_(n), m_(m), puzzle_(n_ * m_) //Build the puzzle with pointers to each square.
        {
            for (size_t i = 0; i < n_; ++i)
                for (size_t j = 0; j < m_; ++j)
                    (*this)(i, j) = s.at(i * m_ + j) - '0';
        }

        u8 & operator () (size_t i, size_t j)  //Overload the & operator to return an array of size i * 3 + j.
        { 
            return puzzle_[i * m_ + j]; 
        }

        u8 const & operator () (size_t i, size_t j) const  //Operator overload to retrieve the value of a tile.
        { 
            return const_cast<Puzzle&>(*this)(i, j); 
        }

        bool operator == (Puzzle const & o) const  //Overload the == operator to compare the equality of two puzzles
        { 
            return puzzle_ == o.puzzle_; 
        
        }

        //Function finds the 0 tile and generates new puzzle states for each direction the tile can move. Returns a vector of new puzzle states.
        vector<Puzzle> Neighbours() const 
        {
            vector<Puzzle> r;
            for (ptrdiff_t i = 0; i < n_; ++i)      //Use a pair of for loops to check the puzzle for the 0 tile
                for (ptrdiff_t j = 0; j < m_; ++j)
                    if ((*this)(i, j) == 0)         //When it finds the 0 tile:
                    {
                        for (pair<int, int> p: vector<pair<int, int>>{{0, -1}, {0, 1}, {-1, 0}, {1, 0}}) //Create a vector of the four possible moves.
                        {
                            ptrdiff_t const ni = i + p.first, nj = j + p.second; //Check that each move results in a valid puzzle state.
                            if (ni < 0 || ni >= n_ || nj < 0 || nj >= m_)
                                continue;
                            Puzzle o = *this;
                            swap(o(i, j), o(ni, nj));  //Make the valid tile move.
                            r.push_back(move(o));      //Push the new puzzle onto vector r.
                        }

                        break;
                    }

            return move(r);     //Return the vector of valid puzzles after moving the 0 tile.
        }

        //Creates a string from the puzzle for use in printing.
        string Str(bool newline = false) const 
        {
            string r;
            for (size_t i = 0; i < n_; ++i)             //Go through each row.
            {
                for (size_t j = 0; j < m_; ++j)         //Go through each column.
                    r.append(1, (*this)(i, j) + '0');   //Apend each item.

                if (newline && i + 1 < n_)              //Every three items, make a new line.
                    r.append(1, '\n');
            }
            return r;                                   //Return the string.
        }

        //Finds the minimum distance between two puzzles.
        size_t MinDist(Puzzle const & to) const 
        {
            size_t r = 0;  //Initialize size(total distance) r to 0.
            for (ptrdiff_t i = 0; i < n_; ++i)      //Iterate through every row and column.
                for (ptrdiff_t j = 0; j < m_; ++j) 
                {
                    auto const v = (*this)(i, j);   //Retrieve a tile and store it in v.

                    if (v == 0)  //If it's the 0 tile, continue to the next tile.
                        continue;

                    size_t dist = size_t(-1);  //Initialize the distance to -1 for error catching.

                    for (ptrdiff_t i2 = 0; i2 < n_; ++i2) 
                    {
                        for (ptrdiff_t j2 = 0; j2 < m_; ++j2)
                            if (to(i2, j2) == v)  //Find the matching tile in the second puzzle and calculate the Manhattan distance.
                            {
                                dist = abs(i - i2) + abs(j - j2);
                                break;
                            }

                        if (dist != size_t(-1))
                            break;
                    }

                    if (dist == -1)
                    {
                        cout << "ERROR\n";
                        return -1;
                    }

                    r += dist;  //Add the tile distance to the total distance.
                }

            return r;  //Return the total distance.
        }

    private:
        size_t n_ = 0, m_ = 0;
        vector<u8> puzzle_;
};

/*
  Parameters: Puzzle const & start - the initial puzzle state
              Puzzle const & goal  - the goal puzzle to be reached
  Returns: The path taken to the goal state.
  Description: Uses the A* algorithm to find the shortest path to the goal state.
*/
vector<Puzzle> AStarSolve(Puzzle const & start, Puzzle const & goal) 
{
    using IdT = string;
    struct Entry 
    {
        Puzzle puzzle;
        size_t gscore = size_t(-1), fscore = size_t(-1);
        IdT came_from{};
    };

    unordered_map<IdT, Entry> entries;  //Stores entries with string types and Entry values.
    map<size_t, set<IdT>> open_set;     //Stores sets of string (IdT) values with a size as key.
    
    //Defines a lambda expression that takes an Entry returns the distance to the goal state.
    auto H = [&](Entry const & e)
    {
        return e.puzzle.MinDist(goal);
    };
    
    {   //Create an entry in data structure "entries" and add it to the open set with puzzle
        //set to start and gscore set to 0. The fscore is then calculated with heuristic H.
        Entry first
        {
            .puzzle = start, .gscore = 0
        };

        first.fscore = H(first);
        entries[first.puzzle.Str()] = first;
        open_set[first.fscore].insert(first.puzzle.Str());
    }
    
    //Function builds a path from a puzzle to the goal state.
    function<vector<Puzzle>(IdT const &, size_t)> ReconstructPath = [&](IdT const & id, size_t depth)
    {
        thread_local vector<Puzzle> path;
        if (id == IdT{})
            return path;

        if (depth == 0)
            path.clear();

        auto const & e = entries.at(id);
        path.insert(path.begin(), e.puzzle);
        return ReconstructPath(e.came_from, depth + 1);
    };
    
    //While the set of vectors is not empty, calculate the best path to the goal state. Return the path.
    while (!open_set.empty()) 
    {
        auto const min_fscore = open_set.begin()->first;
        auto const min_entries = open_set.begin()->second;

        for (auto const & id: min_entries)
            if (entries.at(id).puzzle == goal)
                return ReconstructPath(id, 0);

        open_set.erase(min_fscore);

        for (auto const & cid: min_entries) 
        {
            auto const & cure = entries.at(cid);

            for (auto const & nbid: cure.puzzle.Neighbours()) 
            {
                size_t const tentative_gscore = cure.gscore + 1;

                auto const nid = nbid.Str();
                auto it = entries.find(nid);

                bool is_new = it == entries.end();

                if (is_new || tentative_gscore < it->second.gscore) 
                {
                    if (is_new)
                        it = entries.insert({nid, Entry{.puzzle = nbid}}).first;

                    it->second.came_from = cid;
                    it->second.gscore = tentative_gscore;

                    if (!is_new) 
                    {
                        auto it2 = open_set.find(it->second.fscore);
                        if (it2 != open_set.end() && it2->second.count(nid)) 
                        {
                            it2->second.erase(nid);
                            if (it2->second.empty())
                                open_set.erase(it2);
                        }
                    }

                    it->second.fscore = tentative_gscore + H(it->second);
                    open_set[it->second.fscore].insert(nid);
                }
            }
        }
    }

    cout << "The puzzle can't be solved!\n";  //Display an error if the calculation to find the goal path failed.
    vector<Puzzle> temp;
    return temp;
}

/*
  Parameters: string const & start - the initial puzzle state
              string const & goal  - the goal puzzle to be reached
  Returns: The path taken to the goal state.
  Description: Calls the A* algorithm and prints the result.
*/
int Solve(string const & start, string const & goal) 
{
    auto const v = AStarSolve(start, goal);
    size_t constexpr per_line = 1;
    bool last = false;
    int counter = 0;

    for (size_t i = 0; !last; ++i) 
    {
        for (size_t j = 0; j < n; ++j) 
        {
            for (size_t i2 = 0; i2 < per_line; ++i2) 
            {
                counter++;
                size_t const k = i * per_line + i2;
                if (k >= v.size()) 
                {
                    last = true;
                    for (size_t l = 0; l < (m + 5); ++l)
                        cout << " ";
                } 
                else 
                {
                    auto const & e = v.at(k);
                    auto const s = e.Str(true);
                    size_t pos = 0;

                    for (size_t ip = 0; ip < j; ++ip)
                        pos = s.find('\n', pos) + 1;

                    size_t pos2 = min<size_t>(s.size(), s.find('\n', pos));
                    cout << s.substr(pos, pos2 - pos) << (j == (n / 2) && k + 1 < v.size() ? "    " : "    ");
                }

                cout << (i2 + 1 >= per_line ? "\n" : "");
            }
        }
        cout << endl;
    }
    return counter/4;
}