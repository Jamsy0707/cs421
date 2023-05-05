//James Sieben 200455325
//CS421 A2 8-puzzle problem

#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
using namespace std;

void arrayInput         (int[]);
void arrayPrint         (int[]);
int  checkBestPuzzle    (int, int[], int, int[]);
bool compareArrays      (int[], vector<int[]>&);
void move               (char, int, int[]);
void solve              (int, int[], vector<vector<int>>&);
int  tilesOutOfPlace    (int[]);
int  zeroPosition       (int[]);

const int PUZZLELENGTH = 9;

/*
  Parameters: none
  Returns: 0
  Description: Driver function
*/
int main ()
{
    //char a;
    int count = 1;                                          //Keeps track of number of states. Initialized to one to represent user-entered state.
    vector<int> puzzle;                                     //Integer vcetor holds current puzzle position.
    vector<vector<int>> prevPuzzle;                         //Vector of vectors used to track state path to prevent infinite loops where the next best state is equal to the previous state.
    int puzzleSize = sizeof(puzzle);                        //Represents the size of the puzzle in bytes. Used for copying arrays.

    cout << "This is an implementation of the 8-puzzle problem with the \"Tiles out of place\" heuristic" << endl
         << "Please enter an initial state (Enter 0 for the blank spot):" << endl;

    arrayInput (puzzle);                                    //Call for user to enter values.

    cout << "The sequence is:\n";
    arrayPrint (puzzle);                                    //Print the initial state.
    while (tilesOutOfPlace (puzzle) != 0)                   //Continue to solve until all elements are correct.
    {
        solve (puzzleSize, puzzle, prevPuzzle);             //Call the solve function to find the new puzzle state.
        arrayPrint (puzzle);                                //Print the new puzzle state.
        count++;                                            //Increase the number of states counter.
        if (count > 100)                                    //Prevent the program from running indefinitely in case the puzzle is unsolvable.
            break;
        //cin >> a;
    }

    if (count > 100)
        cout << "The puzzle can't be solved or can't be solved in under 100 moves.\n";
    else
        cout << "The number of states in the sequence is: " << count << "\n\n";

    return 0;
}

/*
  Parameters: int arr[]  - array to be filled
  Returns: nothing
  Description: prompts user for inputs to fill the initial puzzle array
*/
void arrayInput (int arr[])
{
    bool check[PUZZLELENGTH];                               //Bool array is used to check that every value is represented once and only once.
    for (int i = 0; i < PUZZLELENGTH; i++)
        check[i] = false;
    
    for (int i = 0; i < PUZZLELENGTH; i++)                  //Iterate through the entire array.
    {
        cin >> arr[i];                                      //Take in a user input.
        if (arr[i] < 0 || arr[i] >= PUZZLELENGTH)           //If the input is invalid, decrease loop counter, inform the user, and allow them to try again.
        {
            i--;
            cout << "Error: Entry not valid, try again." << endl;
        }
        else if (check[arr[i]] == true)                     //If the input has been used already, decrease loop counter, inform the user, and allow them to try again.
        {
            i--;
            cout << "Error: Values can only be used once." << endl;
        }
        else
            check[arr[i]] = true;                           //Check off the value from the check array.
    }
    cout << endl;
    return;
}

/*
  Parameters: int arr[]  - array to be printed
  Returns: nothing
  Description: prints out all elements of specified array
*/
void arrayPrint (int arr[])
{
    for (int i = 0; i < PUZZLELENGTH; i++)                  //Iterate through the entire array.
    {
        cout << arr[i] << "  ";                             //Two spaces between elements.
        if (i % 3 == 2)                                     //New line every three elements.
            cout << endl;
    }
    cout << endl;
    return;
}

/*
  Parameters: int puzzleSize          - size of puzzle array in bytes
              int tempPuzzle          - copy of puzzle array to be worked on
              int bestTilesOutOfPlace - number of tiles out of place in the current best state
              int bestPuzzle[]        - copy of current best puzzle array state
  Returns: the value of bestTilesOutOfPlace regardless of change
  Description: checks if tempPuzzle is better than bestPuzzle. If it is, make bestPuzzle equal to tempPuzzle.
*/
int checkBestPuzzle (int puzzleSize, int tempPuzzle[], int bestTilesOutOfPlace, int bestPuzzle[])
{
    int tempTilesOutOfPlace = tilesOutOfPlace (tempPuzzle); //Call tilesOutOfPlace to find the number of tiles out of place in the temp state.

    if (tempTilesOutOfPlace < bestTilesOutOfPlace)          //If temp state has less tiles out of place than the current best.
    {
        bestTilesOutOfPlace = tempTilesOutOfPlace;          //Make the temp tiles out of place the new best tiles out of place.
        memcpy(bestPuzzle, tempPuzzle, puzzleSize);         //Make the temp state the new current best state.
    }
    
    return bestTilesOutOfPlace;                             //Return the number of tiles out of place in the best state back to the calling function.
}

/*
  Parameters: int arr1[]        - array to be checked
              vector prevPuzzle - vector of arrays representing previous puzzle states
  Returns: false if arrays are no match is found and true if one is
  Description: compares every element of an array to arrays in a vector to check for equality
*/
bool compareArrays (int arr[], vector<array<int, PUZZLELENGTH>>& prevPuzzle)
{
    /*
    for (int i = 0; i < PUZZLELENGTH; i++)
        if (arr1[i] != arr2[i])
            return false;
    return true;
    */

    for (array<int, PUZZLELENGTH> i : prevPuzzle)   //i is each array in the vector
    {
        int j = 0;                                  //j points to current position in array

        for (auto x : i)                            //auto increment through each position
        {
            if (x != arr[j])                        //if there's a difference, break to next array
                break;
            j++;
            if (j >= 9)                             //If j reaches 9, a full match has been found
                return true;
        }
    }
    return false;                                    //If no match is found, return true
}

/*
  Parameters: char dir     - character representing the direction the zero tile is to be moved
              int zeroPos  - the position of the zero tile in puzzle[]
              int puzzle[] - the current puzzle state
  Returns: nothing
  Description: swaps the zero tile with its neighbour in the specified direction
*/
void move (char dir, int zeroPos, int puzzle[])
{
    int numDir;                                             //Used in the calculation to find where to move the zero tile.

    switch (dir)                                            //Determine which direction to move the tile and how far.
    {
        case 'u': numDir = -3; break;
        case 'd': numDir = 3;  break;
        case 'l': numDir = -1; break;
        case 'r': numDir = 1;  break;
        default: cout << "ERROR in move function." << endl; return;
    }

    puzzle[zeroPos] = puzzle[zeroPos+numDir];               //First, move non-zero tile to zero tile's spot.
    puzzle[zeroPos+numDir] = 0;                             //Second, move zero tile to old non-zero spot.

    return;
}

/*
  Parameters: int puzzleSize    - size of puzzle[] in bytes
              int puzzle[]      - the current puzzle state
              vector prevPuzzle - the previous puzzle states
  Returns: nothing
  Description: tries moving the zero tile in every direction before calling the checkBestPuzzle function to determine the best state to go to next
*/
void solve (int puzzleSize, int puzzle[], vector<array<int, PUZZLELENGTH>>& prevPuzzle)
{
    int  bestTilesOutOfPlace = 999;                         //Used for the tilesOutOfPlace function to determine which state to choose next. Initialized to 999 so first state found is assigned as best to begin.
    int  bestPuzzle[PUZZLELENGTH];                          //Used to store the current best puzzle state when calculating possible next states.
    int  tempPuzzle[PUZZLELENGTH];                          //Stores a copy of the puzzle for use in calculating states without altering the puzzle.
    char charDir;                                           //Stores the direction the zero tile is to move (u, d, l, r).

    int zeroPos = zeroPosition(puzzle);                     //Check where the zero tile is in the current puzzle and save its position.

    for (int i = 0; i < 4; i++)                             //Loop is used so each possible move (u, d, l, r) is tried once.
    {
        switch (i)
        {
            case 0: charDir = 'r'; if ((zeroPos == 2) || (zeroPos == 5) || (zeroPos == 8)) {i++;} else {break;}   //Check if move right is possible. If it is, do it. Else, try the next move.
            case 1: charDir = 'l'; if ((zeroPos == 0) || (zeroPos == 3) || (zeroPos == 6)) {i++;} else {break;}   //Check if move left  is possible. If it is, do it. Else, try the next move.
            case 2: charDir = 'd'; if ((zeroPos == 6) || (zeroPos == 7) || (zeroPos == 8)) {i++;} else {break;}   //Check if move down  is possible. If it is, do it. Else, try the next move.
            case 3: charDir = 'u'; if ((zeroPos == 0) || (zeroPos == 1) || (zeroPos == 2)) {continue;} else {break;}   //Check if move up    is possible. If it is, do it. Else, try the next move.
            default: cout << "Error in switch statement.\n";
        }                          

        memcpy (tempPuzzle, puzzle, puzzleSize);            //Copy puzzle state so it can be worked on without unwanted alterations.
        move (charDir, zeroPos, tempPuzzle);                //Move zero tile in direction specified by previous switch statement.
        if (compareArrays (tempPuzzle, prevPuzzle))         //If the move results in the previous state being achieved (i.e. an infinite loop), don't perform the move.
            continue;
        bestTilesOutOfPlace = checkBestPuzzle (puzzleSize, tempPuzzle, bestTilesOutOfPlace, bestPuzzle);    //Call checkBestPuzzle to see if the tile move resulted in a better state than other moves. Assign the return to keep track of tiles out of place in the current best state.
    }

    prevPuzzle.push_back(puzzle);                           //Copy the current puzzle as the previous puzzle for comparison later.
    memcpy(puzzle, bestPuzzle, puzzleSize);                 //Make the current puzzle state equal to the best puzzle state found by trying every possible move.

    return;
}

/*
  Parameters: int puzzle[]   - the current puzzle state
  Returns: the number of tiles out of place in puzzle[]
  Description: iterates through the puzzle[] array to calulate the number of tiles out of place
*/
int tilesOutOfPlace (int puzzle[])
{
    int x, y, p, q, cost = 0;

    for (int i = 0; i < PUZZLELENGTH; i++)
    {
        switch (i)                                          //Calculate the current coordinates
        {
            case 0: x = 0; y = 0; break;
            case 1: x = 1; y = 0; break;
            case 2: x = 2; y = 0; break;
            case 3: x = 0; y = 1; break;
            case 4: x = 1; y = 1; break;
            case 5: x = 2; y = 1; break;
            case 6: x = 0; y = 2; break;
            case 7: x = 1; y = 2; break;
            case 8: x = 2; y = 2; break;
            default: "Error in x-y coordinate finder in tilesOutOfPlace.\n";
        }

        switch (puzzle[i])                                  //Calculate the final state coordinates
        {
            case 1: p = 0; q = 0; break;
            case 2: p = 1; q = 0; break;
            case 3: p = 2; q = 0; break;
            case 4: p = 0; q = 1; break;
            case 5: p = 1; q = 1; break;
            case 6: p = 2; q = 1; break;
            case 7: p = 0; q = 2; break;
            case 8: p = 1; q = 2; break;
            case 0: p = 2; q = 2; break;
            default: "Error in p-q coordinate finder in tilesOutOfPlace.\n";
        }

        cost += abs(x - p) + abs(y - q);                    //Use the Manhattan distance to calculate cost
    }

    return cost;
}

/*
  Parameters: int puzzle[]   - the current puzzle state
  Returns: the array position of the zero tile
  Description: iterates through the puzzle[] array to find the location of the zero tile
*/
int zeroPosition (int puzzle[])
{
    for (int i = 0; i < PUZZLELENGTH; i++)                  //Iterate through each array position
        if (puzzle[i] == 0)                                 //If element is equal to zero, that position is made the zero position.
            return i;                                       //Return the position of zero.
    return -1;                                         
}