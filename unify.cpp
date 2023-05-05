// James Sieben 200455325
// Unification algorithm implementation

#include <iostream>
#include <string>
#include <cstring>
#include <ctype.h>
using namespace std;


string apply (string E, string SUBS1)                      //Function "apply" applies S2 to S1 and returns the result as a string
{
    if (SUBS1.find("/") == string::npos)                   //If S2 is not in the form A/B return the original string
        return {E};

    string toReplace = SUBS1.substr(0, SUBS1.find("/"));   //The portion of the string before '/' is to be replaced
    string replacement = SUBS1.substr(SUBS1.find("/")+1);  //The portion of the string after '/' is replacing

    size_t pos = 0;
    while ((pos = E.find(toReplace, pos)) != string::npos) //Go through the string and find all instances to replace
    {
        E.replace(pos, toReplace.length(), replacement);   //Replace 
        pos += replacement.length();                       //Skip to after the portion that has been replaced to avoid needless checking 
    }

    return {E};                                            //Return modified string
}

string unify (string E1, string E2)                        //Function "unify" implements the unification algo
{
    bool isFunction = false;
    int i;

    if ((E1[0] == '(') || (E1[0] == ')'))                  //If first term is a parenthese we get an infinite loop so simply remove it
        E1.erase(0,1);

    if ((E2[0] == '(')  || (E2[0] == ')'))                 //Same premise with the second string
        E2.erase(0,1);

    for (i = 0; i < E1.length(); i++)                      //If string E1 contains an invalid character, return "Invalid"
        if (!(isalpha(E1[i]) || (E1[i] == '/') || (E1[i] == '(') || (E1[i] == ')') || (E1[i] == ',')))
            return {"Invalid"};

    for (i = 0; i < E2.length(); i++)                      //If string E2 contains an invalid character, return "Invalid"
        if (!(isalpha(E2[i]) || (E2[i] == '/') || (E2[i] == '(') || (E2[i] == ')') || (E2[i] == ',')))
            return {"Invalid"};

    for (i = 0; i < E1.length(); i++)                      //If the string has brackets, mark it as a function to skip the next steps
        if ((E1[i] == '(') || (E1[i] == ')'))
        {
            isFunction = true;
            break;
        }
    if (isFunction == false)
        for (i = 0; i < E2.length(); i++)                  //If the string has brackets, mark it as a function to skip the next steps
            if ((E2[i] == '(') || (E2[i] == ')'))
            {
                isFunction = true;
                break;
            }

    if ((isFunction == false) &&                           //If the string starts with a lower case character, return an empty string if they're equal or "FAIL" if they're not
        (E1.empty() || (E1[0] >= 'a' && E1[0] <= 'z')) &&
        (E2.empty() || (E2[0] >= 'a' && E2[0] <= 'z')))
    {
        if (E1 == E2)
            return {};
        else
            return {"FAIL"};
    }

    if (E1[0] >= 'A' && E1[0] <= 'Z')                      //If the first string starts with an upper case character, return "FAIL" if E1 occurs in E2 or E2/E1 if not
    {
        if (E2.find(E1) != string::npos)
            return {"FAIL"};
        else
            return {E2+"/"+E1};
    }

    else if (E2[0] >= 'A' && E2[0] <= 'Z')                //If the seond string starts with an upper case character, return "FAIL" if E2 occurs in E1 or E1/E2 if not
    {
        if (E1.find(E2) != string::npos)
            return {"FAIL"};
        else
            return {E1+"/"+E2};
    }

    else if (E1.empty() || E2.empty())                   //If either string is empty, return "FAIL"
        return {"FAIL"};

    else
    {
        string HE1(1, E1[0]);                            //HE1 = the first character of E1
        string HE2(1, E2[0]);                            //HE2 = the first character of E2

        string SUBS1 = unify(HE1, HE2);                  //Unify the heads of E1 and E2 and return FAIL if the unification of those heads returns FAIL
        if (SUBS1 == "FAIL")
            return {"FAIL"};

        string TE1 = apply(E1.substr(1), SUBS1);         //Apply the unification of the heads to the rest of the strings E1 and E2
        string TE2 = apply(E2.substr(1), SUBS1);

        string SUBS2 = unify(TE1, TE2);                  //Unify the results of the last step and return "FAIL" if "FAIL" is returned from the unification
        if (SUBS2 == "FAIL")
            return {"FAIL"};
        else                                             //Else apply SUBS2 to SUBS1
        {
            return {apply(SUBS1, SUBS2)};
        }
    }

    return {"FAIL"};                                     //Returns "FAIL" if the end of the function is reached without returning anything else
}

int main()                                               //The function "main" takes two strings from the user and prints out the unification of the two
{
    string E1, E2;
    cout << "This is an implementation of the unification algorithm.\n";
    while ((E1 != "q") && (E2 != "q"))                   //If the user enters 'q', quit the program
    {
        cout << "Please enter the first term: ";
        cin >> E1;
        cout << "Please enter the second term: ";
        cin >> E2;

        string result = unify(E1, E2);

        if (result == "FAIL")                            //If "FAIL" has been returned, no unification exists
            cout << "There is no unifier.\n\n";
        else if (result == "Invalid")                    //If "Invalid" has been returned, there was a problem with the strings entered
            cout << "Invalid.\n\n";
        else
        {
            string first = result.substr(result.find("/")+1);   //Look for everything before the '/' in the returned string
            string second = result.substr(0, result.find("/")); //Look for everything after the '/' in the returned string

            cout << "The unifier is: " << first << " = " << second << "\n\n";
        }
    }

    return 0;
}