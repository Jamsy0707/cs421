/*
    James Sieben 200455325
    CS 421, Assignment 3
    Bayesian Network Solver Using Variable Elimination

    To Use:
        Compile and run as normal. The program is usually slow to start/calculate due to its length.
        When entering a query, do not type "Pr(A|B)," instead type "A|B" without quotes.
        There is an issue with the code, calculating values changes the values for future calculations,
           to avoid this problem the program must be rerun for each query.

    Works Cited:
        This code was adapted from the GitHub project titled "bayonet" by user mpatacchiola.
*/


#include <algorithm>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <initializer_list>
#include <list>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>


namespace bayonet{

    class JointProbabilityTable {

        public:
            JointProbabilityTable();
            JointProbabilityTable(std::vector<unsigned int> variablesTotStatesVector);
            ~JointProbabilityTable();

            std::pair<std::vector<unsigned int>, double> ReturnRow(unsigned int index);

            double ReturnMarginal(unsigned int variableIndex, unsigned int variableState);
            std::vector<unsigned int> ReturnKey(unsigned int index);
            double GetProbability(std::vector<unsigned int> variablesStatesVector);
            bool SetProbability(std::vector<unsigned int> variablesStatesVector, double probability);
            bool AddToProbability(std::vector<unsigned int> variablesStatesVector, double valueToAdd);
            void Print();
            void PrintProbability(std::vector<unsigned int> variablesStatesVector);
            void PrintMarginals();
            void PrintMarginal(unsigned int variableIndex, int count, std::string num);

            void NormalizeProbabilities(double alpha=0);
            void RandomizeProbabilities();
            void ResetProbabilities(double valueToSet=0);

            void AddVariable(unsigned int totStates);
            void Clear();

            unsigned int ReturnRowsNumber();
            const std::map<std::vector<unsigned int>,double>& ReturnJointMap();


        private:
            std::vector<unsigned int> mVariablesTotStatesVector;
            std::map<std::vector<unsigned int>,double> mJointMap;

            void FillMap(std::vector<unsigned int> variablesTotStatesVector);

    };

    class MarginalProbabilityTable {

        public:
            MarginalProbabilityTable(unsigned int totVariables, unsigned int totStates);
            MarginalProbabilityTable(std::vector<unsigned int> variablesTotStatesVector);
            ~MarginalProbabilityTable();

            bool SetProbability(unsigned int variableIndex, unsigned int stateIndex, double probability);
            bool AddToProbability(unsigned int variableIndex, unsigned int stateIndex, double probability);
            double GetProbability(unsigned int variableIndex, unsigned int stateIndex);

            unsigned int ReturnMostProbableState(unsigned int variableIndex);

            bool SetProbabilities(unsigned int index, std::vector<double> probabilitiesVector);
            std::vector<double> GetProbabilities(unsigned int index);

            void ResetProbabilities(double valueToSet=0);
            void NormalizeProbabilities();

            void Print();
            void PrintVariable(unsigned int index);

        private:
            std::vector<std::vector<double>> marginalTable;

    };

    class ConditionalProbabilityTable {

        public:
            ConditionalProbabilityTable(unsigned int NodeStatesNumber);
            ConditionalProbabilityTable(unsigned int NodeStatesNumber, std::vector<unsigned int> parentsStatesList);
            ~ConditionalProbabilityTable();

            std::pair<std::vector<unsigned int>, std::vector<double>> ReturnRow(unsigned int index);
            std::vector<unsigned int> ReturnParentsState(unsigned int index);

            std::vector<unsigned int> FindParentState(unsigned int parentIndex, unsigned int parentState);

            double GetProbability(unsigned int variableState, std::vector<unsigned int> parentsStates);
            std::vector<double> GetProbabilities(std::vector<unsigned int> parentsStates);
            std::vector<double> GetProbabilities(unsigned int variableState, unsigned int parentIndex, unsigned int parentsState);
            bool SetProbabilities(std::vector<unsigned int> parentsStates, std::vector<double> probabilities);
            bool AddToProbability(unsigned int variableState, std::vector<unsigned int> parentsStates, double valueToAdd);
            void Print();
            void PrintProbabilities(std::vector<unsigned int> parentsStates);

            void NormalizeProbabilities();
            void RandomizeProbabilities();

            unsigned int ReturnSample(std::vector<unsigned int> parentsStates);

            void ResetProbabilities(double valueToSet=0.0);
            void AddVariable(unsigned int totStates);
            void Clear();

            unsigned int ReturnRowsNumber();
            unsigned int ReturnColumnsNumber();

        private:
            std::vector<unsigned int> mTotalParentsStates;
            std::map<std::vector<unsigned int>,std::vector<double> > conditionalMap;
            void FillMap(unsigned int NodeStatesNumber, std::vector<unsigned int> parentsStates);

    };

    class Bayesnode {

        public:

            /// Enum class colour
            enum colour {
            WHITE, ///< is used for unobserved nodes
            GREY, ///< is is used for partially observed nodes
            BLACK, ///< is used for observed nodes
            };

            ConditionalProbabilityTable conditionalTable;

            Bayesnode(unsigned int numberOfStates);
            ~Bayesnode();

            int ReturnNumberOfStates();

            void SetLabel(std::string);
            std::string GetLabel();
            void SetNumericLabel(int);
            int GetNumericLabel();

            void SetColour(colour);
            colour GetColour();

            bool AddToAdjacencyList(unsigned int);
            bool RemoveFromAdjacencyList(unsigned int);
            const std::list<unsigned int>&  ReturnAdjacencyList();
            bool IsInAdjacencyList(unsigned int);
            unsigned int SizeOfAdjacencyList();

            bool SetEvidence(int evidenceState);
            unsigned int GetEvidence();
            bool IsEvidence();


        private:
            int mEvidence;
            colour mCurrentColour;
            int mNumberOfStates;
            std::string mNodeLabel;
            int mNodeNumericLabel;
            std::list<unsigned int> adjacencyList;

    };

    class Bayesnet{

        public:

            Bayesnet(unsigned int totNodes, unsigned int totStates);
            Bayesnet(std::vector<unsigned int> nodesTotStatesVector);
            ~Bayesnet();

            Bayesnode& operator[](unsigned int index);

            bool AddEdge(unsigned int FirstNode, unsigned int SecondNode);
            bool RemoveEdge(unsigned int FirstNode, unsigned int SecondNode);
            bool HasEdge(unsigned int FirstNode, unsigned int SecondNode);

            unsigned int ReturnNumberOfNodes();
            unsigned int ReturnNumberOfEdges();

            std::list<unsigned int> ReturnOutEdges(unsigned int index);
            std::list<unsigned int> ReturnInEdges(unsigned int index);
            unsigned int ReturnNumberOutEdges(unsigned int index);
            unsigned int ReturnNumberInEdges(unsigned int index);

            std::list<unsigned int> ReturnTopologicalList();
            std::list<unsigned int> ReturnRootList();
            std::list<unsigned int> ReturnLeafList();

            std::vector<unsigned int> ReturnTotalStates();
            std::vector<unsigned int> ReturnNotEvidenceNodes();
            std::vector<unsigned int> ReturnEvidenceNodes();

            double GetNodeProbability(unsigned int index, std::vector<unsigned int> variablesStatesVector);

            void ResetAllColours();

            bool IsTree();
            //bool IsPolytree(); //TODO
            bool IsMultiConnected();

            bool IsRoot(unsigned int);
            bool IsLeaf(unsigned int); //TODO
            //unsigned int ReturnMarkovBlanketSize(unsigned int index); //TODO
            //double ReturnAverageMarkovBlanketSize(unsigned int index); //TODO

            const std::vector<Bayesnode>& ReturnNodesVector();
            void FillJointProbabilityTable();

            std::list<unsigned int> BreadthFirstSearch(unsigned int startingNode);
            std::list<unsigned int> DepthFirstSearch(unsigned int startingNode);

        private:
            std::vector<Bayesnode> nodesVector;
            std::pair<std::list<unsigned int>, unsigned int> RawDepthFirstSearch(unsigned int index);

    };

    class MaximumLikelihoodLearning {

        public:
            MaximumLikelihoodLearning();
            ~MaximumLikelihoodLearning();

            Bayesnet ReturnUpdatedNetwork(Bayesnet net, std::vector<std::vector<unsigned int>>& trainingDataset);

        private:

    };

    class GibbsSampler {

        public:
            GibbsSampler();
            ~GibbsSampler();

            std::vector<unsigned int> ReturnSample(bayonet::Bayesnet& net);
            std::vector<unsigned int> ReturnSample(bayonet::Bayesnet& net, std::vector<unsigned int> startingVector);

            std::vector<std::vector<unsigned int>> AccumulateSamples(Bayesnet& net, unsigned int cycles);

            void PrintSample(bayonet::Bayesnet& net, unsigned int cycles = 1);
            JointProbabilityTable ReturnJointProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles);
            MarginalProbabilityTable ReturnMarginalProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles);

        private:
            int ReturnInteger(const int minRange, const int maxRange);

    };

    class LWSampler {

        public:
            LWSampler();
            ~LWSampler();

            std::pair< std::vector<unsigned int>, double> ReturnSample(bayonet::Bayesnet& net);

            std::vector<std::pair<std::vector<unsigned int>, double>> AccumulateSamples(Bayesnet&, unsigned int cycles);

            void PrintSample(bayonet::Bayesnet& net, unsigned int cycles = 1);
            JointProbabilityTable ReturnJointProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles);
            MarginalProbabilityTable ReturnMarginalProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles);


        private:

    };

    class RejectionSampler {

        public:
            RejectionSampler();
            ~RejectionSampler();
            std::vector<unsigned int> ReturnSample(bayonet::Bayesnet& net);
            std::vector<std::vector<unsigned int>> AccumulateSamples(Bayesnet&, unsigned int cycles);
            std::vector<std::vector<unsigned int>> AccumulateAndDiscardSamples(Bayesnet&, unsigned int cycles);
            void PrintSample(bayonet::Bayesnet& net, unsigned int cycles = 1);
            JointProbabilityTable ReturnJointProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles);
            MarginalProbabilityTable ReturnMarginalProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles);

        private:

    };

    Bayesnet::Bayesnet(unsigned int totNodes, unsigned int totStates){
        nodesVector.reserve(totNodes);
        for(unsigned int i=0; i<totNodes; i++){
            Bayesnode my_node(totStates);
            nodesVector.push_back(my_node); //filling the nodes vector
        }
    }

    Bayesnet::Bayesnet(std::vector<unsigned int> nodesTotStatesVector) {
        nodesVector.reserve(nodesTotStatesVector.size());
        for(auto it=nodesTotStatesVector.begin(); it!=nodesTotStatesVector.end(); ++it){
            Bayesnode my_node(*it);
            nodesVector.push_back(my_node); //filling the nodes vector
        }
    }

    /**
    * It destroys the object.
    *
    **/
    Bayesnet::~Bayesnet(){}

    /**
    * Operator overload [] it is used to return thereference to the node stored inside the net
    * It is possible to access the methods of the single node in a easier way.
    * Example: net[2].IsRoot();  // It checks if the third node is a root node
    * @param index the number of the element stored inside the net
    * @return it returns a reference to the node
    **/
    Bayesnode& Bayesnet::operator[](unsigned int index){
        if (index >= nodesVector.size()) throw std::domain_error("Error: Out of Range index.");
        return nodesVector[index];
    }

    /**
    * It adds an Edge between two nodes.
    *
    * @param firstNode the parent node
    * @param secondNode the child node
    **/
    bool Bayesnet::AddEdge(unsigned int firstNode, unsigned int secondNode){
        if(firstNode > nodesVector.size() || secondNode > nodesVector.size()){
            std::cerr << "ERROR: Out of range index" << std::endl;
            return false;
        }
        if(firstNode == secondNode) return false;
        unsigned int first_tot_states = nodesVector[firstNode].ReturnNumberOfStates();
        nodesVector[firstNode].AddToAdjacencyList(secondNode);
        nodesVector[secondNode].conditionalTable.AddVariable(first_tot_states);
        return true;
    }

    /**
    * It removes an Edge between two nodes.
    *
    * @param firstNode the parent node
    * @param secondNode the child node
    **/
    bool Bayesnet::RemoveEdge(unsigned int firstNode, unsigned int secondNode){
        if(firstNode > nodesVector.size() || secondNode > nodesVector.size()){
            std::cerr << "ERROR: Out of range index" << std::endl;
            return false;
        }
        if(firstNode == secondNode) return false;
        nodesVector[firstNode].RemoveFromAdjacencyList(secondNode);
        return true;
    }

    /**
    * It checks if an Edge between two nodes exist.
    *
    * @param firstNode the parent node
    * @param secondNode the child node
    **/
    bool Bayesnet::HasEdge(unsigned int firstNode, unsigned int secondNode){
    if(firstNode > nodesVector.size() || secondNode > nodesVector.size()){
    std::cerr << "ERROR: Out of range index" << std::endl;
    return false;
    }
    return nodesVector[firstNode].IsInAdjacencyList(secondNode);
    }

    /**
    * It returns the number of nodes.
    *
    **/
    unsigned int Bayesnet::ReturnNumberOfNodes(){
    return nodesVector.size();
    }

    /**
    * It returns the number of edges.
    *
    **/
    unsigned int Bayesnet::ReturnNumberOfEdges(){
    unsigned int tot_edges = 0;
    for(auto it=nodesVector.begin(); it!=nodesVector.end(); ++it){
    tot_edges += it->SizeOfAdjacencyList();
    }

    return tot_edges;
    }

    /**
    * It returns output index list
    *
    **/
    std::list<unsigned int> Bayesnet::ReturnOutEdges(unsigned int index){
    std::list<unsigned int> temp_list;
    temp_list = nodesVector[index].ReturnAdjacencyList();
    return temp_list;
    }

    /**
    * It returns the input index list.
    *
    **/
    std::list<unsigned int> Bayesnet::ReturnInEdges(unsigned int index){
    std::list<unsigned int> temp_list;
    unsigned int counter = 0;
    for(auto it=nodesVector.begin(); it!=nodesVector.end(); ++it){
    if((*it).IsInAdjacencyList(index) == true) temp_list.push_back(counter);
    counter++;
    }
    return temp_list;
    }

    /**
    * It returns the number of Incoming edges from the node specified in the index.
    *
    **/
    unsigned int Bayesnet::ReturnNumberOutEdges(unsigned int index){
    return ReturnOutEdges(index).size();
    }

    /**
    * It returns the number of ingoing edges from the node specified in the index.
    *
    **/
    unsigned int Bayesnet::ReturnNumberInEdges(unsigned int index){
    return ReturnInEdges(index).size();
    }


    /**
    * The topological sort algorithm creates a linear ordering of the vertices. If edge (u,v) appears in the graph, then u comes before v in the ordering. 
    * The topological order is used during the sampling phase for sending queries to a node only when all the parents were already sorted.
    * The algorithm for topological ordering use the DepthFirstSearch function to calculate the finish time of every node.
    * Each node in a DAG goes from a node of higher finish time to a node of lower node finish time. 
    * The problem with cyclic graph is due to the fact that back edges go from nodes of lower finish time to nodes of higher finish time 
    *
    * @return It returns a list of index to nodes sorted in topological order.
    **/
    std::list<unsigned int> Bayesnet::ReturnTopologicalList() {

    std::list<unsigned int> list_to_return;
    std::multimap<unsigned int, unsigned int> index_time_map;

    //Iterating through each node and applying the DFS algorithm to calculate the finish time
    for(unsigned int nodes_counter = 0; nodes_counter < nodesVector.size(); nodes_counter++){
    //std::shared_ptr<std::list<unsigned int>> spToList;
    //auto deep_list = std::make_shared<std::list<unsigned int>>();
    auto deep_list = DepthFirstSearch(nodes_counter);
    index_time_map.insert(std::pair<unsigned int,unsigned int>(deep_list.size(), nodes_counter));
    }

    //filling the output list from the multimap
    for(auto it_map=index_time_map.begin(); it_map!=index_time_map.end(); ++it_map){
    list_to_return.push_front(it_map->second);
    }

    return list_to_return;
    }

    /**
    * I return a list of integers, containing the index of
    * all the Root nodes.
    *
    **/
    std::list<unsigned int> Bayesnet::ReturnRootList(){
    std::list<unsigned int> list_to_return;
    for(unsigned int i=0; i<nodesVector.size(); i++){
    if(IsRoot(i) == true) list_to_return.push_back(i);
    }
    return list_to_return;
    }

    /**
    * I return a list of integers, containing the index of
    * all the Leaf nodes.
    *
    **/
    std::list<unsigned int> Bayesnet::ReturnLeafList(){
    std::list<unsigned int> list_to_return;
    for(unsigned int i=0; i<nodesVector.size(); i++){
    if(IsLeaf(i) == true) list_to_return.push_back(i);
    }
    return list_to_return;
    }

    /**
    * Given the the index of a variable and the states of all
    * the variables, it returns the associated probability.
    *
    * @param variableState
    * @param parentsStates
    **/
    double Bayesnet::GetNodeProbability(unsigned int index, std::vector<unsigned int> variablesStatesVector){

    auto in_list = ReturnInEdges(index);
    std::vector<unsigned int> key_vector;

    for(auto it_list=in_list.begin(); it_list!=in_list.end(); ++it_list){
    key_vector.push_back(variablesStatesVector[*it_list]);
    }

    unsigned int variable_state = variablesStatesVector[index];
    double result = nodesVector[index].conditionalTable.GetProbability(variable_state, key_vector);
    return result;
    }

    /**
    * Reset all the nodes colours to white.
    * 
    **/
    void Bayesnet::ResetAllColours(){
    for(auto it_node=nodesVector.begin(); it_node!=nodesVector.end(); ++it_node){
    it_node->SetColour(Bayesnode::colour::WHITE);
    }
    }

    /**
    * It returns a list with the total number of states for each node.
    * 
    * @return it returns the list of total states
    **/
    std::vector<unsigned int> Bayesnet::ReturnTotalStates(){

    std::vector<unsigned int> vector_to_return;
    for(auto it_node=nodesVector.begin(); it_node!=nodesVector.end(); ++it_node){
    vector_to_return.push_back( it_node->ReturnNumberOfStates() );
    }
    return vector_to_return;
    }

    /**
    * It returns a list containing the index to the nodes which are not Evidence nodes.
    * 
    * @return it returns the list of not evidence node
    **/
    std::vector<unsigned int> Bayesnet::ReturnNotEvidenceNodes(){
    std::vector<unsigned int> list_to_return;
    unsigned int counter = 0;
    for(auto it_node=nodesVector.begin(); it_node!=nodesVector.end(); ++it_node){
    if(it_node->IsEvidence() == false) list_to_return.push_back( counter );
    counter++;
    }
    return list_to_return;
    }

    /**
    * It returns a list containing the index to the nodes which are Evidence nodes.
    * 
    * @return it returns the list of evidence node
    **/
    std::vector<unsigned int> Bayesnet::ReturnEvidenceNodes(){
    std::vector<unsigned int> list_to_return;
    unsigned int counter = 0;
    for(auto it_node=nodesVector.begin(); it_node!=nodesVector.end(); ++it_node){
    if(it_node->IsEvidence() == true) list_to_return.push_back( counter );
    counter++;
    }
    return list_to_return;
    }

    /**
    * A Bayesian network is a Tree if it has only
    * one root node and no cycles.
    * 
    * @return it returns true if the network is a tree
    **/
    bool Bayesnet::IsTree(){
    if(ReturnRootList().size() == 1 && IsMultiConnected()==false) return true;
    else return false;
    }

    /**
    * Multi-Connected Bayesian network is a network with
    * multiple possible paths from a starting node to an
    * arriving node. To check if the network is multi-
    * connected, a DepthFirst search is done from the
    * root nodes, during the DFS traversal it is possible
    * to see if some nodes were already visited.
    * 
    * @return it returns true if the network is multi-conncted
    **/
    bool Bayesnet::IsMultiConnected(){
    ResetAllColours();
    unsigned int tot_neurons = nodesVector.size();
    for(unsigned int i_node=0; i_node<tot_neurons; i_node++){ 
    if(IsRoot(i_node) == true){
    auto deep_pair = RawDepthFirstSearch(i_node);
    if (deep_pair.second > 0) return true;
    }
    }
    return false;
    }

    /**
    * A root node is a node without parents but with children.
    * 
    * @return it returns true if the node is a root node.
    **/
    bool Bayesnet::IsRoot(unsigned int index){
    auto out_list = ReturnOutEdges(index);
    auto in_list = ReturnInEdges(index);
    if(in_list.size() == 0 && out_list.size() > 0) return true;
    else return false;
    }

    /**
    * A leaf node is a node without children but with parents.
    * 
    * @return it returns true if the node is a leaf node.
    **/
    bool Bayesnet::IsLeaf(unsigned int index){
    auto out_list = ReturnOutEdges(index);
    auto in_list = ReturnInEdges(index);
    if(in_list.size() > 0 && out_list.size() == 0) return true;
    else return false;
    }

    /**
    * It return a const reference to the internal vector used to store the nodes.
    * 
    * @return it returns a const reference
    **/
    const std::vector<Bayesnode>& Bayesnet::ReturnNodesVector(){
    return nodesVector;
    }


    /**
    * Breadth First Search algorithm. Starting from a node it performs a breadth-first traversal of the network.
    *
    * @param startingNode the index of the node 
    * @return it returns a list of index representing the order of nodes visited
    **/
    std::list<unsigned int> Bayesnet::BreadthFirstSearch(unsigned int startingNode){

    std::list<unsigned int> list_to_return;

    //Mark all nodes as white (not visited)
    ResetAllColours();

    // Create a queue for BFS
    std::list<int> queue;

    // Mark the current node as visited and enqueue it
    nodesVector[startingNode].SetColour(Bayesnode::colour::BLACK);
    queue.push_back(startingNode);

    while(!queue.empty())
    {
    unsigned int dequeue_node = queue.front();
    //std::cout << dequeue_node << " ";
    list_to_return.push_back(dequeue_node);
    queue.pop_front();

    // Get all adjacent vertices of the dequeued node.
    // If a adjacent has not been visited, then mark it visited and enqueue it.
    auto temp_list = ReturnOutEdges(dequeue_node);

    for(auto it = temp_list.begin(); it != temp_list.end(); ++it)
    {

    if(nodesVector[*it].GetColour() == Bayesnode::colour::WHITE){
        nodesVector[*it].SetColour(Bayesnode::colour::BLACK);
        queue.push_back(*it);
    }
    }
    }
    return list_to_return;
    }


    /**
    * Depth First Search algorithm. Starting from a node it performs a depth-first traversal of the network.
    *
    * @param startingNode the index of the node
    * @param spToList a shared pointer to a list, that is filled recursively by the algorithm.
    * @param resetColour if true it reset all the colour to white before starting the algorithm
    **/
    std::list<unsigned int> Bayesnet::DepthFirstSearch(unsigned int startingNode) {
    ResetAllColours();
    auto deep_pair = RawDepthFirstSearch(startingNode);
    auto deep_list = deep_pair.first;
    deep_list.reverse();
    return deep_list;
    }


    /**
    * Depth First Search algorithm. It is a raw version that returns also
    * the number of nodes that have a multiple connection.
    * It is used also from IsMultiConnected()
    *
    * @param startingNode the index of the node
    **/
    std::pair<std::list<unsigned int>, unsigned int> Bayesnet::RawDepthFirstSearch(unsigned int startingNode){

    nodesVector[startingNode].SetColour(Bayesnode::colour::GREY);
    std::list<unsigned int> list_to_return;
    auto temp_list = ReturnOutEdges(startingNode);
    unsigned int multiconnections = 0;

    //Chain iteration to all the children
    for(auto it = temp_list.begin(); it != temp_list.end(); ++it)
    {
    if(nodesVector[*it].GetColour() == Bayesnode::colour::BLACK){
        //std::cout << "  MULTICONNECTED HERE: " << *it;
        multiconnections += 1;
    }
    if(nodesVector[*it].GetColour() == Bayesnode::colour::WHITE){
        nodesVector[*it].SetColour(Bayesnode::colour::GREY);
        auto children_pair = RawDepthFirstSearch(*it);
        multiconnections += children_pair.second;   
        list_to_return.merge(children_pair.first);
    }
    }

    // Done Visiting StartingNode
    nodesVector[startingNode].SetColour(Bayesnode::colour::BLACK);
    list_to_return.push_back(startingNode);
    return std::make_pair(list_to_return, multiconnections);
    }

    Bayesnode::Bayesnode(unsigned int numberOfStates) : conditionalTable(numberOfStates), mEvidence(-1), mCurrentColour(WHITE), mNumberOfStates(numberOfStates), mNodeLabel(""), mNodeNumericLabel(0)
    {
    }

    Bayesnode::~Bayesnode()
    {
    }

    /**
    * It returns the number of states associated with the node.
    * 
    * @return the number of states associated with the node
    **/
    int Bayesnode::ReturnNumberOfStates(){
    return mNumberOfStates;
    }

    /**
    * It sets the label associated with the node.
    * 
    * @param nodeLabel the string to associate with the label
    **/
    void Bayesnode::SetLabel(std::string nodeLabel){
    mNodeLabel = nodeLabel;
    }

    /**
    * It returns the label associated with the node.
    * 
    * @return the label associated with the node
    **/
    std::string Bayesnode::GetLabel(){
    return mNodeLabel;
    }

    /**
    * It sets the numeric label of the node.
    * 
    * @return the label associated with the node
    **/
    void Bayesnode::SetNumericLabel(int numericLabel){
    mNodeNumericLabel = numericLabel;
    }

    /**
    * It returns the numeric label associated with the node.
    * 
    * @return the numeric label associated with the node
    **/
    int Bayesnode::GetNumericLabel(){
    return mNodeNumericLabel;
    }

    /**
    * It sets the colour of the node.
    * The colour is used in searching algorithms.
    * White vertices are undiscovered. 
    * Grey is an adjacent vertex that is not already discovered.
    * Black vertices are discovered and are adjacent to only other black or gray vertices.
    * 
    **/
    void Bayesnode::SetColour(colour colourToSet){
    mCurrentColour = colourToSet;
    }

    /**
    * It return the colour associated with the node.
    * 
    **/
    Bayesnode::colour Bayesnode::GetColour(){
    return mCurrentColour;
    }

    /**
    * Adding an index to the adjacency list of the node.
    * 
    * @param index
    **/
    bool Bayesnode::AddToAdjacencyList(unsigned int index){

    //looking for copies
    for(std::list<unsigned int>::iterator it = adjacencyList.begin(); it != adjacencyList.end(); ++it) {
    if( *it == index) return false;
    }
    adjacencyList.push_back(index);
    //conditionalTable.AddVariable(totStates);

    return true;
    }


    /**
    * Removing an index to the adjacency list of the node.
    * 
    * @param index
    * @return it returns true if the element was found and erased
    **/
    bool Bayesnode::RemoveFromAdjacencyList(unsigned int index){
    bool result = false;
    //looking if present
    for(std::list<unsigned int>::iterator it = adjacencyList.begin(); it != adjacencyList.end(); ++it) {
    if( *it == index) result = true;;
    }
    adjacencyList.remove(index);
    return result;
    }


    /**
    * It returns a const reference to the internal adjacency list.
    *
    * @return a const reference to the list
    * 
    **/
    const std::list<unsigned int>& Bayesnode::ReturnAdjacencyList(){
    return adjacencyList;
    }

    /**
    * It looks for the index specified in input inside the adjacency list.
    *
    * @param index 
    * @return true if the element is inside the list
    * 
    **/
    bool Bayesnode::IsInAdjacencyList(unsigned int index){
    //looking for edges
    for(std::list<unsigned int>::iterator it = adjacencyList.begin(); it != adjacencyList.end(); ++it) {
    if( *it == index) return  true;
    }
    return false;
    }

    /**
    * It returns the number of Incoming edges
    *
    * @return 
    * 
    **/
    unsigned int Bayesnode::SizeOfAdjacencyList(){
    return adjacencyList.size();
    }

    /**
    * A node is an evidence when is outcome state is given.
    * 
    * @return it returns -1 if the node is not an evidence.
    * Otherwise it return the state sets as evidence.
    **/
    bool Bayesnode::IsEvidence(){
    if(mEvidence < 0) return false;
    else return true;
    }

    /**
    * A node is an evidence when is outcome state is given.
    *
    * @param evidenceState the state to set as evidence.
    * @return it returns true if the state was correctly set.
    * 
    **/
    bool Bayesnode::SetEvidence(int evidenceState){
    if(evidenceState < ReturnNumberOfStates()){
    mEvidence = evidenceState;
    return true;
    }else{
    return false;
    }
    }

    /**
    * It returns the evidence. If the node is not an evidence node it returns -1
    *
    * @return it returns true if the state was correctly set.
    * 
    **/
    unsigned int Bayesnode::GetEvidence(){
    return mEvidence;
    }

    ConditionalProbabilityTable::ConditionalProbabilityTable(unsigned int NodeStatesNumber){

    //Wrong input check
    if(NodeStatesNumber < 2) NodeStatesNumber=2;

    //Fill the map
    FillMap(NodeStatesNumber, {});

    }

    /**
    * The constructor of the object.
    *
    * @param NodeStatesNumber is the number of states of the node associated with the table
    * @param parentsStates this is a list that represents the number of states of each parent of the node
    * If the node has 3 boolean states as parents the list is {2, 2, 2}
    **/
    ConditionalProbabilityTable::ConditionalProbabilityTable(unsigned int NodeStatesNumber, std::vector<unsigned int> parentsStates){

    //Wrong input check
    if(NodeStatesNumber < 2) NodeStatesNumber=2;

    //Fill the map
    FillMap(NodeStatesNumber, parentsStates);

    //Set the Parents state vector
    mTotalParentsStates = parentsStates;
    }

    /**
    * Destroying the object
    *
    **/
    ConditionalProbabilityTable::~ConditionalProbabilityTable(){}

    /**
    * Return a single row of the table.
    * The row is a pair containing the parent state
    * and the probabilities for all the variable states.
    *
    * @param index
    **/
    std::pair<std::vector<unsigned int>, std::vector<double>> ConditionalProbabilityTable::ReturnRow(unsigned int index){
    std::vector<unsigned int> parent_vector;
    std::vector<double> state_vector;
    if(index > conditionalMap.size()) return make_pair(parent_vector, state_vector); //out of range > empty vector returned
    auto it_map=conditionalMap.begin();
    std::advance(it_map, index);
    parent_vector = it_map->first;
    state_vector = it_map->second;
    return make_pair(parent_vector, state_vector);
    }

    /**
    * Return the Parents state at a given index
    *
    * @param index
    **/
    std::vector<unsigned int> ConditionalProbabilityTable::ReturnParentsState(unsigned int index){
    std::vector<unsigned int> vector_to_return;
    if(index > conditionalMap.size()) return vector_to_return; //out of range > empty vector returned
    auto it_map=conditionalMap.begin();
    std::advance(it_map, index);
    vector_to_return = it_map->first;
    return vector_to_return;
    }

    /**
    * It find in which rows the specified parent has the specified state.
    * It returns a list of integers, representing the positions of the
    * parent in the CPT.
    *
    * @param parentIndex
    * @param parentState
    **/
    std::vector<unsigned int> ConditionalProbabilityTable::FindParentState(unsigned int parentIndex, unsigned int parentState){
    unsigned int row = 0;
    std::vector<unsigned int> rows_vector;

    //Iterate through map content:
    for (auto it_map=conditionalMap.begin(); it_map!=conditionalMap.end(); ++it_map){
    std::vector<unsigned int> parents_vector = it_map->first;

    if(parents_vector.size()==0){
    std::cerr << "ERROR: the Conditional Table does not have any parent" << std::endl;
    return rows_vector;
    }
    if(parentIndex > parents_vector.size()-1){
    std::cerr << "ERROR: parentIndex out of range in Conditional Table" << std::endl;
    return rows_vector;
    }
    //std::cout  << "  "  << "    CHEK 1.1.1.1.CON " << parentIndex << " # " << parents_vector.size() << " " << row << std::endl;
    if(parents_vector[parentIndex] == parentState) rows_vector.push_back(row);

    row++;
    }
    return rows_vector;
    }

    /**
    * Given the state of the variable and a vector key it returns the associated probability
    *
    * @param variableState
    * @param parentsStates
    **/
    double ConditionalProbabilityTable::GetProbability(unsigned int variableState, std::vector<unsigned int> parentsStates){
    auto row_vector = conditionalMap[parentsStates];
    //std::cout << "AT: " << row_vector.at(variableState) << std::endl;
    return row_vector.at(variableState);
    }



    /**
    * Given a vector key it returns the associated probabilities
    *
    * @param parentsStates
    **/
    std::vector<double> ConditionalProbabilityTable::GetProbabilities(std::vector<unsigned int> parentsStates){
    return conditionalMap[parentsStates];
    }

    /**
    * Given a variable State, a specific parent index and a parent state,
    * it returns all the probabilities associated with that configuration
    * finding them through the rows of the CPT
    *
    * @param variableState
    * @param parentIndex
    * @param parentsState
    **/
    std::vector<double> ConditionalProbabilityTable::GetProbabilities(unsigned int variableState, unsigned int parentIndex, unsigned int parentsState){
    std::vector<double> vector_to_return;
    //Iterate through map content:
    for (auto it_map=conditionalMap.begin(); it_map!=conditionalMap.end(); ++it_map){
    if(it_map->first[parentIndex] == parentsState) vector_to_return.push_back( it_map->second[variableState] );
    }
    return vector_to_return;
    }

    /**
    * Given a vector key it sets the probabilities associated
    *
    * @param parentsStates
    * @param probabilities
    **/
    bool ConditionalProbabilityTable::SetProbabilities(std::vector<unsigned int> parentsStates, std::vector<double> probabilities){
    if(conditionalMap.find(parentsStates) != conditionalMap.end()){
    conditionalMap[parentsStates] = probabilities;
    return true;
    }else{
    return false;
    }
    }

    /**
    * Given a vector key and a variable state, it adds a value
    * to the associated probability.
    *
    * @param variableState
    * @param parentsStates
    * @param valueToAdd
    **/
    bool ConditionalProbabilityTable::AddToProbability(unsigned int variableState, std::vector<unsigned int> parentsStates, double valueToAdd){
    if(variableState > conditionalMap.at(parentsStates).size()){
    std::cerr << "ERROR: Conditional Table out of range index." << std::endl;
    return false; 
    }

    conditionalMap.at(parentsStates)[variableState] += valueToAdd;
    return true;
    }

    /**
    * It prints the conditional table on the terminal.
    * If the number of columns and rows is huge, the terminal could cut parts of the output.
    *
    **/
    void ConditionalProbabilityTable::Print(){
    
    std::cout << std::endl;

    //std::cout << '|' << std::setw(5) << "KEY" << std::setw(5) << '|' << std::endl;

    //Iterate through map content:
    for (auto it_map=conditionalMap.begin(); it_map!=conditionalMap.end(); ++it_map){

    //std::cout << "KEY: ";
    std::cout << '|' << std::setw(5);
    for (auto it_key=it_map->first.begin(); it_key!=it_map->first.end(); ++it_key){
        std::cout << (*it_key);
        if(it_key!=it_map->first.end()-1) std::cout << "-";
    }
    std::cout << std::setw(5);

    //std::cout << "  PROB: ";
    for (auto it_data=it_map->second.begin(); it_data!=it_map->second.end(); ++it_data){
        std::cout << "| "  << std::setw(6)  << (*it_data) << std::setw(6);  //<< std::setprecision(4)
    }

    std::cout << '|' << std::endl;
    }

    std::cout << std::endl << std::endl;
    std::cout << "COLUMNS .... " << ConditionalProbabilityTable::ReturnColumnsNumber() << std::endl;
    std::cout << "ROWS    .... " << ConditionalProbabilityTable::ReturnRowsNumber() << std::endl;
    }

    /**
    * It prints the probabilities associated with a particular set of parents.
    * If the number of columns is huge, the terminal could cut parts of the output.
    *
    * @param parentsStates
    **/
    void ConditionalProbabilityTable::PrintProbabilities(std::vector<unsigned int> parentsStates){

    std::vector<double> prob_vector = conditionalMap[parentsStates];

    std::cout << '\n';

    //std::cout << "KEY: ";
    std::cout << '|' << std::setw(5);
    for (auto it=parentsStates.begin(); it!=parentsStates.end(); ++it){
    std::cout << *it;
    if(it!=parentsStates.end()-1) std::cout << "-";
    }
    std::cout << std::setw(5);

    //std::cout << "  PROB: ";
    for (auto it=prob_vector.begin(); it!=prob_vector.end(); ++it){
    std::cout << "| "  << std::setw(6) << *it << std::setw(6);
    }

    std::cout << '|' << '\n';
    }

    /**
    * It normalize the probabilities inside each row of the table.
    *
    **/
    void ConditionalProbabilityTable::NormalizeProbabilities(){
    
    //Iterate through map content:
    for (auto it_map=conditionalMap.begin(); it_map!=conditionalMap.end(); ++it_map){
    //define an accumulator variable
    //and assigne to it the value of zero
    double accumulator = 0;
    //Iterate through the vector of probabilities for accumulate the sum
    for (auto it_data=it_map->second.begin(); it_data!=it_map->second.end(); ++it_data){
        accumulator += *it_data;
    }

    //Iterate through the vector of probabilities for the normalization
    for (auto it_data=it_map->second.begin(); it_data!=it_map->second.end(); ++it_data){
        //If the accumulator is zero, then all the variables are zero
        //It prevent a problem in case of a zero divisor
        if(accumulator != 0) *it_data = *it_data / accumulator;
        else *it_data = 0.0;
    }
    }
    }

    /**
    * It randomize the probabilities inside each row of the table.
    * The probabilities are also normalized.
    **/
    void ConditionalProbabilityTable::RandomizeProbabilities(){

    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_real_distribution<> real_dist(0, 1);


    //Iterate through map content:
    for (auto it_map=conditionalMap.begin(); it_map!=conditionalMap.end(); ++it_map){

    //save the size of the vector and then clear it
    unsigned int numberOfStates = it_map->second.size();
    it_map->second.clear();
    double accumulator = 0;

    //generate random values and push them inside the state vector
    for (unsigned int i = 0; i < numberOfStates; i++) {
    double temp_value = real_dist(generator);
    it_map->second.push_back(temp_value);
    accumulator += temp_value;      
    }

    //Iterate through the vector of probabilities for the normalization
    for (auto it_data=it_map->second.begin(); it_data!=it_map->second.end(); ++it_data){
    *it_data = *it_data / accumulator;
    }
    }

    }

    /**
    * It returns a random sample, given a particular configuration of the parents
    *
    * @param parentsStates
    **/
    unsigned int ConditionalProbabilityTable::ReturnSample(std::vector<unsigned int> parentsStates){
    std::random_device random_device;
    std::mt19937 generator(random_device());

    std::discrete_distribution<unsigned int> values_distribution (conditionalMap[parentsStates].begin(), conditionalMap[parentsStates].end());

    return values_distribution(generator);
    }

    /**
    * It add a new parent to the table.
    *
    * @param totStates the total number of states to assign to the new parent
    **/
    void ConditionalProbabilityTable::AddVariable(unsigned int totStates=2){

    //Check for wrong input
    if(totStates <2) totStates=2;

    //Push the new total in the vector
    mTotalParentsStates.push_back(totStates);

    unsigned int tot_columns = ReturnColumnsNumber();

    //clear the map
    conditionalMap.clear();

    FillMap(tot_columns, mTotalParentsStates);
    }

    /**
    * It set to zero all the entries of the table.
    *
    * @param valueToSet
    **/
    void ConditionalProbabilityTable::ResetProbabilities(double valueToSet){
    //Iterate through map content:
    for (auto it_map=conditionalMap.begin(); it_map!=conditionalMap.end(); ++it_map){

    //Iterate through the vector of probabilities
    for (auto it_data=it_map->second.begin(); it_data!=it_map->second.end(); ++it_data){
    *it_data = valueToSet;
    }
    }
    }

    /**
    * It clear the content of the table.
    *
    **/
    void ConditionalProbabilityTable::Clear(){
    conditionalMap.clear();
    mTotalParentsStates.clear();
    }

    /**
    * It returns the total number of rows.
    *
    **/
    unsigned int ConditionalProbabilityTable::ReturnRowsNumber(){
    return conditionalMap.size();
    }

    /**
    * It returns the total number of columns.
    * The number of columns coincides with the number of states.
    *
    **/
    unsigned int ConditionalProbabilityTable::ReturnColumnsNumber(){
    if(conditionalMap.size() == 0) return 0;
    else{
    return conditionalMap.begin()->second.size();
    }
    }

    /**
    * Private function. It's a low level function for filling the content of the map.
    * It works like an odometer and create all the possible combination of states starting from the parents states.
    *
    * @param NodeStatesNumber
    * @param parentsStates
    **/
    void ConditionalProbabilityTable::FillMap(unsigned int NodeStatesNumber, std::vector<unsigned int> parentsStates){
    //0- In case of a root node without parents, it adds only one row
    if(parentsStates.size() == 0){
    std::vector<double> data_vector;
    std::vector<unsigned int> key_vector;
    for(unsigned int i=0; i<NodeStatesNumber; i++){
    data_vector.push_back((double)1.0 / (double)NodeStatesNumber);
    }
    conditionalMap.insert( std::make_pair( key_vector, data_vector ) );
    return;
    }

    //1- Creation of the vector that contains the Subsets of unsigned int
    std::vector<std::vector<unsigned int>> vector_container; 
    for (auto it = parentsStates.begin() ; it != parentsStates.end() ; ++it)
    {
    std::vector<unsigned int> temp_vector;
    unsigned int total_states = (*it);
    for(unsigned int state_counter=0; state_counter<total_states; state_counter++){
    temp_vector.push_back(state_counter);
    }
    vector_container.push_back(temp_vector);
    }

    //2- Creation of the vector that contains the Iterators
    std::vector< std::vector<unsigned int>::iterator > iterator_container;
    for (auto it = vector_container.begin() ; it != vector_container.end() ; ++it){
    std::vector<unsigned int>::iterator temp_iterator;
    temp_iterator = (*it).begin();
    iterator_container.push_back(temp_iterator);
    }

    //3- filling the data_vector
    std::vector<double> data_vector;
    data_vector.reserve(NodeStatesNumber);
    for(unsigned int p_counter=0; p_counter<NodeStatesNumber; p_counter++){
    data_vector.push_back((double) 1.0 / (double)NodeStatesNumber);
    }

    //4- Cascade iteration for storing the whole set of combination
    unsigned int K = iterator_container.size();

    while (iterator_container[0] != vector_container[0].end()) {

        //Clear the key_vector before pushing a new key
        std::vector<unsigned int> key_vector;
        key_vector.reserve(parentsStates.size());

        //Filling the key_vector 
        for (auto it_key=iterator_container.begin(); it_key!=iterator_container.end(); ++it_key){    
        key_vector.push_back( (*(*it_key)) ); //pushing a new key
        }

        //Assign a new row in the conditional map
        conditionalMap.insert( std::make_pair( key_vector, data_vector ) );

    //Increment the odometer by 1 
    //ATTENTION: It must called at the end of the while otherwise the first key get missed...
    ++iterator_container[K-1]; 
    for (int i = K-1; (i > 0) && (iterator_container[i] == vector_container[i].end()); --i) {
        //subtracting the counter
        iterator_container[i] = vector_container[i].begin();
        ++iterator_container[i-1];
    }
    }
    }

    GibbsSampler::GibbsSampler(){}

    GibbsSampler::~GibbsSampler(){}

    /**
    * It returns a single sample picking up it from the Bayesian network
    *
    * @param net the Bayesian network to use for picking up the sample.
    *
    **/
    std::vector<unsigned int> GibbsSampler::ReturnSample(bayonet::Bayesnet& net){

    //Declaring the vector to return
    std::vector<unsigned int> vector_to_return;

    //Declaring the local variables
    auto topo_list = net.ReturnTopologicalList();
    std::map<unsigned int, std::pair<bool, unsigned int>> sample_map;
    unsigned int tot_nodes = net.ReturnNumberOfNodes();

    //Fill the sample_map and the vector_to_return with zeros
    for(unsigned int i=0; i<tot_nodes; i++){
    auto my_pair = std::make_pair<bool, unsigned int>(false,0);
    auto map_pair = std::make_pair(i, my_pair);
    sample_map.insert(map_pair);
    vector_to_return.push_back(0);
    }

    //Cycle through the topological list
    for(auto it_topo=topo_list.begin(); it_topo!=topo_list.end(); ++it_topo){
    auto in_list = net.ReturnInEdges(*it_topo);
    std::vector<unsigned int> key_vector;
    //Cycle through the in list for creating the key
    for(auto it_in=in_list.begin(); it_in!=in_list.end(); ++it_in){
    auto map_value = sample_map[*it_in]; //return the pair stored inside the map
    bool value_check = map_value.first; //return the boolean of the pair
    unsigned int value_sample = map_value.second; //return the unsigned int of the pair
    if(value_check == false){
        std::cerr << "EXCEPTION: the topological order was not respected!" << std::endl;   
    }else{
        key_vector.push_back(value_sample); //push   
    }
    }
    //Key completed, asking for the sample
    auto sp_to_node = net[*it_topo];
    //Checking if the node is Evidence
    if(sp_to_node.IsEvidence() == true){
    //unsigned int node_evidence = sp_to_node->GetEvidence();
    std::pair<bool,unsigned int> pair_to_store = std::make_pair<bool, unsigned int>(true, sp_to_node.GetEvidence());
    sample_map[*it_topo] = pair_to_store;
    }else{
    //If it is not an Evidence then storing the sample in the local index vector
    std::pair<bool,unsigned int> pair_to_store = std::make_pair<bool, unsigned int>(true, sp_to_node.conditionalTable.ReturnSample(key_vector));
    sample_map[*it_topo] = pair_to_store;
    }

    }

    //Fill the vector_to_return with the result stored inside the map
    //This cycle restore the nodes order (not the topological order)
    for(unsigned int i=0; i<tot_nodes; i++){
    auto temp_pair = sample_map[i];
    vector_to_return[i] = temp_pair.second; //Storing the sample into the vector to return
    }

    return vector_to_return;
    }

    /**
    * It returns a single sample picking up it from the Bayesian network
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param startingVector it is the vector with the current states of all the nodes
    * it is used as starting point for the Markov chain sampler
    **/
    std::vector<unsigned int> GibbsSampler::ReturnSample(bayonet::Bayesnet& net, std::vector<unsigned int> startingVector){
    std::vector<unsigned int> vector_to_return;

    //Return list of not-evidence nodes
    auto not_evidence_nodes_list = net.ReturnNotEvidenceNodes();

    //Checking for empty Networks
    if(net.ReturnNumberOfNodes() == 0) return vector_to_return; //The net is empty

    //Check if the network contains only evidence nodes
    if(not_evidence_nodes_list.size() == 0){
    return startingVector;
    }

    unsigned int temp_int = ReturnInteger(0, not_evidence_nodes_list.size()-1);
    unsigned int selected_node = not_evidence_nodes_list[temp_int]; //ReturnRandomFromList(not_evidence_nodes);
    auto sp_to_node = net[selected_node];

    //It is necessary to iterate the different states of the node
    //in order to build the final conditional distribution
    std::vector<unsigned int> temp_sample_vector = startingVector;
    unsigned int node_tot_states = sp_to_node.ReturnNumberOfStates();
    std::vector<double> conditional_distribution_vector;
    for(unsigned int i_states=0; i_states<node_tot_states; i_states++){

        temp_sample_vector[selected_node] = i_states;
        //1- The first phase consists in finding the probability
        //of the selected node given its parents
        auto in_list = net.ReturnInEdges(selected_node);
        std::vector<unsigned int> key_vector;
        //Cycle through the in list for creating the key
        for(auto it_in=in_list.begin(); it_in!=in_list.end(); ++it_in){
        unsigned int selected_state = startingVector[*it_in]; //take the state of the parents of the node
            key_vector.push_back(selected_state); //push the state in the key vector  
        }
        //Key completed, asking for the probability
        double parents_prob = sp_to_node.conditionalTable.GetProbability(i_states ,key_vector);

        //2- The second phase consists in finding the probability
        //for each children of the current node, given the parents
        auto out_list = net.ReturnOutEdges(selected_node);
        std::vector<double> children_prob_vector;
        for(auto it_out=out_list.begin(); it_out!=out_list.end(); ++it_out){ 
        double returned_prob = net.GetNodeProbability(*it_out, temp_sample_vector);
        children_prob_vector.push_back(returned_prob);
        }

        //3- The last phase is the multiplication of the different
        //probabilities in order to get the final Markov blanket conditional probability
        double conditional_prob = parents_prob;
        for(auto it_children=children_prob_vector.begin(); it_children!=children_prob_vector.end(); ++it_children){ 
        conditional_prob *= *it_children;
        }
        conditional_distribution_vector.push_back(conditional_prob);
    }

    //4- It is necessary to normalize the conditional distribution before sampling
    double accumulator = 0.0;
    for(auto it_dist=conditional_distribution_vector.begin(); it_dist!=conditional_distribution_vector.end(); ++it_dist){
    accumulator += *it_dist;
    }
    for(auto it_dist=conditional_distribution_vector.begin(); it_dist!=conditional_distribution_vector.end(); ++it_dist){
    *it_dist = *it_dist / accumulator;
    }

    //5- It is time to sample from the distribution
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::discrete_distribution<unsigned int> values_distribution (conditional_distribution_vector.begin(), conditional_distribution_vector.end());
    unsigned int node_sample =values_distribution(generator);
    startingVector[selected_node] = node_sample;

    //6- Returning the modified sample vector
    return startingVector;
    }


    /*
    std::vector<std::vector<unsigned int>> GibbsSampler::AccumulateSamples(Bayesnet& net, unsigned int cycles){
    std::vector<std::vector<unsigned int>> vector_to_return;
    //The first sample is taken at random
    auto starting_sample_vector = ReturnSample(net);
    vector_to_return.push_back(starting_sample_vector);
    //Return list of not-evidence nodes
    auto not_evidence_nodes_list = net.ReturnNotEvidenceNodes();
    //Checking for empty Networks
    if(net.ReturnNumberOfNodes() == 0) return vector_to_return; //The net is empty
    //Check if the network contains only evidence nodes
    if(not_evidence_nodes_list.size() == 0){
    for(unsigned int i=0; i<cycles-1; i++){
    vector_to_return.push_back(starting_sample_vector);
    }
    return vector_to_return;
    }
    //The starting_sample is modified one variable at time
    //and the result is stored. The vairable to change is
    //choosen at random as well.
    for(unsigned int i=0; i<cycles-1; i++){
    unsigned int temp_int = ReturnInteger(0, not_evidence_nodes_list.size()-1);
    unsigned int selected_node = not_evidence_nodes_list[temp_int]; //ReturnRandomFromList(not_evidence_nodes);
    auto sp_to_node = net[selected_node];
    //It is necessary to iterate the different states of the node
    //in order to build the final conditional distribution
    std::vector<unsigned int> temp_sample_vector = starting_sample_vector;
    unsigned int node_tot_states = sp_to_node->ReturnNumberOfStates();
    std::vector<double> conditional_distribution_vector;
    for(unsigned int i_states=0; i_states<node_tot_states; i_states++){
        temp_sample_vector[selected_node] = i_states;
        //1- The first phase consists in finding the probability
        //of the selected node given its parents
        auto in_list = net.ReturnInEdges(selected_node);
        std::vector<unsigned int> key_vector;
        //Cycle through the in list for creating the key
        for(auto it_in=in_list.begin(); it_in!=in_list.end(); ++it_in){
        unsigned int selected_state = starting_sample_vector[*it_in]; //take the state of the parents of the node
            key_vector.push_back(selected_state); //push the state in the key vector  
        }
        //Key completed, asking for the sample
        double parents_prob = sp_to_node->conditionalTable.GetProbability(i_states ,key_vector);
        //Storing the sample in the local index vector
        //unsigned int the_sample = sp_to_node->conditionalTable.ReturnSample(key_vector);
        //starting_sample_vector[selected_node] = the_sample;
        //vector_to_return.push_back(starting_sample_vector);
        //2- The second phase consists in finding the probability
        //for each children of the current node, given the parents
        auto out_list = net.ReturnOutEdges(selected_node);
        std::vector<double> children_prob_vector;
        for(auto it_out=out_list.begin(); it_out!=out_list.end(); ++it_out){ 
        double returned_prob = net.GetNodeProbability(*it_out, temp_sample_vector);
        children_prob_vector.push_back(returned_prob);
        }
        //3- The last phase is the multiplication of the different
        //probabilities in order to get the final Markov blanket conditional probability
        double conditional_prob = parents_prob;
        for(auto it_children=children_prob_vector.begin(); it_children!=children_prob_vector.end(); ++it_children){ 
        conditional_prob *= *it_children;
        }
        conditional_distribution_vector.push_back(conditional_prob);
    }
    //4- It is necessary to normalize the conditional distribution before sampling
    double accumulator = 0.0;
    for(auto it_dist=conditional_distribution_vector.begin(); it_dist!=conditional_distribution_vector.end(); ++it_dist){
    accumulator += *it_dist;
    }
    for(auto it_dist=conditional_distribution_vector.begin(); it_dist!=conditional_distribution_vector.end(); ++it_dist){
    *it_dist = *it_dist / accumulator;
    }
    //5- It is time to sample from the distribution
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::discrete_distribution<unsigned int> values_distribution (conditional_distribution_vector.begin(), conditional_distribution_vector.end());
    unsigned int node_sample =values_distribution(generator);
    starting_sample_vector[selected_node] = node_sample;
    //6- Pushing inside the output vector the new sample
    vector_to_return.push_back(starting_sample_vector);
    }
    return vector_to_return;
    }
    */

    /**
    * This method is different from the same methods in the other samplers.
    * The first sample is obtained at random. The next samples are choosen
    * for each node picking up a value from the Markov blanket of the node.
    * This probability is proportional to the probability of the variable
    * given its parents times the probability of each child given its 
    * respective parents. 
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    *
    **/
    std::vector<std::vector<unsigned int>> GibbsSampler::AccumulateSamples(Bayesnet& net, unsigned int cycles){
    std::vector<std::vector<unsigned int>> accumulated_samples_vector; 
    std::vector<unsigned int> starting_sample_vector = ReturnSample(net);
    
    for(unsigned int i=0; i<cycles; i++){
    accumulated_samples_vector.push_back(starting_sample_vector);
    starting_sample_vector =  ReturnSample(net, starting_sample_vector);
    }
    return accumulated_samples_vector;
    }

    /**
    * It prints the result of the sampling.
    * It is possible to do it for different iterations.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    *
    **/
    void GibbsSampler::PrintSample(bayonet::Bayesnet& net, unsigned int cycles){
    std::vector<unsigned int> starting_sample_vector = ReturnSample(net);
    
    for(unsigned int i=0; i<cycles; i++){
    std::cout << i+1 << " ..... ";
    for(auto it=starting_sample_vector.begin(); it!=starting_sample_vector.end(); ++it){
    std::cout << *it << " ";
    }
    std::cout << std::endl;
    starting_sample_vector =  ReturnSample(net, starting_sample_vector);
    }
    }

    /**
    * It creates a Joint Probability table starting from the Bayesian network and sampling for
    * the number of iterations specified.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    * @return it return a Joint Probability Table object
    **/
    JointProbabilityTable GibbsSampler::ReturnJointProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles){

    //0-Declare the JPT
    JointProbabilityTable joint_table(net.ReturnTotalStates());

    //1-reset JPT
    joint_table.ResetProbabilities();

    //2-Accumulate samples
    auto samples_vector = AccumulateSamples(net, cycles);

    //3-Add sample to JPT
    for(auto it_sample=samples_vector.begin(); it_sample!=samples_vector.end(); ++it_sample){
    joint_table.AddToProbability(*it_sample, 1);
    }

    //4-Normalize JPT
    joint_table.NormalizeProbabilities();

    //5-Return JPT
    return joint_table;
    }

    /**
    * It creates a Marginal Probability table starting from the Bayesian network and sampling for
    * the number of iterations specified.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    * @return it return a Marginal Probability Table object
    **/
    MarginalProbabilityTable GibbsSampler::ReturnMarginalProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles){
    MarginalProbabilityTable marginal_table(net.ReturnTotalStates());

    marginal_table.ResetProbabilities();

    //2-Accumulate samples
    auto samples_vector = AccumulateSamples(net, cycles);

    //2-Accumulate samples and Add sample and weight to JPT
    for(auto it_acc=samples_vector.begin(); it_acc!=samples_vector.end(); ++it_acc){
    auto sample_pair = *it_acc;
    unsigned int var_counter = 0;
    //Iteration through each element of the sample
    //each element is a variable and the value is the state
    for(auto it_sample=sample_pair.begin(); it_sample!=sample_pair.end(); ++it_sample){
    marginal_table.AddToProbability(var_counter, *it_sample, 1.0);
    var_counter++;
    } 

    }

    marginal_table.NormalizeProbabilities();

    return marginal_table;


    }

    /**
    * Return a random generated integer.
    * It is used a uniform distribution.
    * @param minRange it is the minimum value to use for selecting the integer
    * @param maxRange it is the maximum value to use for selecting the integer
    * @return It returns the selected integer
    */
    int GibbsSampler::ReturnInteger(const int minRange, const int maxRange) {
        std::srand (std::time(NULL));
        std::random_device generator_device;
        std::mt19937_64 gen(generator_device());
        std::uniform_int_distribution<int> distribution_device(minRange,maxRange);
        return distribution_device(gen);
    }

    JointProbabilityTable::JointProbabilityTable(){

    }

    /**
    * The constructor of the object.
    *
    * @param variablesTotStates this is a list that represents the number of states of each parent of the node
    * If the node has 3 boolean states as parents the list is {2, 2, 2}
    **/
    JointProbabilityTable::JointProbabilityTable(std::vector<unsigned int> variablesTotStates){

    //Fill the map
    FillMap(variablesTotStates);

    //Set the Parents state vector
    mVariablesTotStatesVector = variablesTotStates;
    }

    /**
    * Destroying the object
    *
    **/
    JointProbabilityTable::~JointProbabilityTable(){}

    /**
    * It returns the specified row from the table
    *
    * @param index
    **/
    std::pair<std::vector<unsigned int>, double> JointProbabilityTable::ReturnRow(unsigned int index){
    std::vector<unsigned int> parent_vector;
    double probability;
    if(index > mJointMap.size()){
    std::cerr << "ERROR: Joint Table out of range index" << std::endl;
    return make_pair(parent_vector, probability); //out of range > empty vector returned
    }
    auto it_map=mJointMap.begin();
    std::advance(it_map, index);
    parent_vector = it_map->first;
    probability = it_map->second;
    return make_pair(parent_vector, probability);
    }

    /**
    * It returns the marginal probability of a given variable.
    *
    * @param variableIndex
    * @param variableState
    **/
    double JointProbabilityTable::ReturnMarginal(unsigned int variableIndex, unsigned int variableState){

        double total = 0.0000;
        //Iterate through map content:
        for (auto it_map=mJointMap.begin(); it_map!=mJointMap.end(); ++it_map){
            if(it_map->first.at(variableIndex) == variableState) total += it_map->second;
        }
        return total;
    }

    /**
    * It returns the key assciated with the index.
    *
    * @param index
    **/
    std::vector<unsigned int> JointProbabilityTable::ReturnKey(unsigned int index){
    unsigned int counter = 0;
    for(auto it=mJointMap.begin(); it!=mJointMap.end(); ++it){
    if(counter == index) return (*it).first;
    counter ++;
    }
    return {};
    }

    /**
    * Given a vector key it returns the associated probabilities
    *
    * @param variablesTotStates
    **/
    double JointProbabilityTable::GetProbability(std::vector<unsigned int> variablesTotStates){
    return mJointMap[variablesTotStates];
    }

    /**
    * Given a vector key it sets the probabilities associated
    *
    * @param variablesTotStates
    * @param probabilities
    **/
    bool JointProbabilityTable::SetProbability(std::vector<unsigned int> variablesTotStates, double probability){
    if(mJointMap.find(variablesTotStates) != mJointMap.end()){
    mJointMap[variablesTotStates] = probability;
    //std::cout << "P: " << probability << std::endl;
    return true;
    }else{
    return false;
    }
    }

    bool JointProbabilityTable::AddToProbability(std::vector<unsigned int> variablesStatesVector, double valueToAdd){
    if(mJointMap.find(variablesStatesVector) != mJointMap.end()){
    mJointMap[variablesStatesVector] += valueToAdd;
    return true;
    }else{
    return false;
    }
    }

    /**
    * It prints the conditional table on the terminal.
    * If the number of columns and rows is huge, the terminal could cut parts of the output.
    *
    **/
    void JointProbabilityTable::Print(){
    
    std::cout << std::endl;

    //std::cout << '|' << std::setw(5) << "KEY" << std::setw(5) << '|' << std::endl;

    //Iterate through map content:
    for (auto it_map=mJointMap.begin(); it_map!=mJointMap.end(); ++it_map){

    //THE KEY
    for (auto it_key=it_map->first.begin(); it_key!=it_map->first.end(); ++it_key){
        std::cout << (*it_key);
        if(it_key!=it_map->first.end()-1) std::cout << "-";
    }

    //THE DATA
    std::cout << " ..... " << it_map->second  << std::endl;
    }

    std::cout << std::endl;
    std::cout << "ROWS ..... " << JointProbabilityTable::ReturnRowsNumber() << std::endl;
    std::cout << std::endl;
    }

    /**
    * It prints the probabilities associated with a particular set of parents.
    * If the number of columns is huge, the terminal could cut parts of the output.
    *
    **/
    void JointProbabilityTable::PrintProbability(std::vector<unsigned int> variablesTotStates){

    //THE KEY
    for (auto it=variablesTotStates.begin(); it!=variablesTotStates.end(); ++it){
    std::cout << *it;
    if(it!=variablesTotStates.end()-1) std::cout << "-";
    }

    //THE DATA
    std::cout << " ..... " << mJointMap[variablesTotStates]  << std::endl;
    std::cout << std::endl;
    }

    /**
    * It prints the marginal probabilities associated with all the variables.
    * If the number of columns is huge, the terminal could cut parts of the output.
    *
    **/
    void JointProbabilityTable::PrintMarginals(){
    unsigned int variable_counter = 0;
    for (auto it=mVariablesTotStatesVector.begin(); it!=mVariablesTotStatesVector.end(); ++it){
    for(unsigned int i=0; i<*it; i++){
    double probability = ReturnMarginal(variable_counter, i);
    std::cout << "Variable: " << variable_counter << " State: " << i << " ..... " << probability << std::endl;
    }
    variable_counter++;
    } 
    }

    /**
    * It prints the marginal probabilities associated with a specific variable.
    * If the number of columns is huge, the terminal could cut parts of the output.
    *
    **/
    void JointProbabilityTable::PrintMarginal(unsigned int variableIndex, int count, std::string num){
        std::cout << std::fixed;
        std::cout << std::setprecision(4);
        unsigned int variableTotStates = mVariablesTotStatesVector[variableIndex];
        if (count == 1){
            for(unsigned int i=0; i<variableTotStates; i++){
                double probability = ReturnMarginal(variableIndex, i);
                std::cout << i << "     " << probability << std::endl;
            }
        }
        else if (count == 2){
            for(unsigned int i=0; i<variableTotStates; i++){
                double probability = ReturnMarginal(variableIndex, i);
                std::cout << num << "     " << i << "     " << probability << std::endl;
            }
        }
        else if (count == 3){
            for(unsigned int i=0; i<variableTotStates; i++){
                double probability = ReturnMarginal(variableIndex, i);
                std::cout << num[0] << "     " << num[1] << "     " << i << "     " << probability << std::endl;
            }
        }
    }

    /**
    * It normalize the probabilities inside each row of the table.
    *
    **/
    void JointProbabilityTable::NormalizeProbabilities(double alpha){
    //double accumulator = 0;

    //If the normalizing constant is less than zero it is not possible to divide
    //then it will be find thorugh iteration
    if(alpha <= 0){
    alpha = 0;
    //Iterate through map content:
    for (auto it_map=mJointMap.begin(); it_map!=mJointMap.end(); ++it_map){
        alpha += it_map->second;
    }
    }

    //Iterate through map content for the normalization
    for (auto it_map=mJointMap.begin(); it_map!=mJointMap.end(); ++it_map){
        it_map->second = it_map->second / alpha;
    }
    }

    /**
    * It randomize the probabilities inside each row of the table.
    * The probabilities are also normalized.
    **/
    void JointProbabilityTable::RandomizeProbabilities(){

    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_real_distribution<> real_dist(0, 1);
    double accumulator = 0;

    //Iterate through map for adding the values
    for (auto it_map=mJointMap.begin(); it_map!=mJointMap.end(); ++it_map){

    //generate random values and push them inside the map
    double temp_value = real_dist(generator);
    it_map->second = temp_value;
    accumulator += temp_value;      
    }

    //Iterate through map for the normalization
    for (auto it_map=mJointMap.begin(); it_map!=mJointMap.end(); ++it_map){
    it_map->second = it_map->second / accumulator;
    }
    }

    /**
    * It reset the probabilities inside each row of the table.
    * Each value is set to zero.
    **/
    void JointProbabilityTable::ResetProbabilities(double valueToSet){
    //Iterate through map for the normalization
    for (auto it_map=mJointMap.begin(); it_map!=mJointMap.end(); ++it_map){
    it_map->second = valueToSet;
    }
    }

    /**
    * It add a new variable to the table.
    *
    * @param totStates the total number of states to assign to the new variable
    **/
    void JointProbabilityTable::AddVariable(unsigned int totStates=2){

    //Check for wrong input
    if(totStates < 2) totStates=2;

    //Push the new total in the vector
    mVariablesTotStatesVector.push_back(totStates);

    //clear the map
    mJointMap.clear();

    FillMap( mVariablesTotStatesVector);
    }

    /**
    * It clear the content of the table.
    *
    **/
    void JointProbabilityTable::Clear(){
    mJointMap.clear();
    mVariablesTotStatesVector.clear();
    }

    /**
    * It returns the total number of rows.
    *
    **/
    unsigned int JointProbabilityTable::ReturnRowsNumber(){
    return mJointMap.size();
    }


    /**
    * Private function. It's a low level function for filling the content of the map.
    * It works like an odometer and create all the possible combination of states starting from the parents states.
    *
    **/
    void JointProbabilityTable::FillMap(std::vector<unsigned int> variablesTotStates){
    //0- In case of a root node without parents, it adds only one row
    if(variablesTotStates.size() == 0){
    return;
    }

    //1- Creation of the vector that contains the Subsets of unsigned int
    std::vector<std::vector<unsigned int>> vector_container; 
    for (auto it = variablesTotStates.begin() ; it != variablesTotStates.end() ; ++it)
    {
    std::vector<unsigned int> temp_vector;
    unsigned int total_states = (*it);
    for(unsigned int state_counter=0; state_counter<total_states; state_counter++){
    temp_vector.push_back(state_counter);
    }
    vector_container.push_back(temp_vector);
    }

    //2- Creation of the vector that contains the Iterators
    std::vector< std::vector<unsigned int>::iterator > iterator_container;
    for (auto it = vector_container.begin() ; it != vector_container.end() ; ++it){
    std::vector<unsigned int>::iterator temp_iterator;
    temp_iterator = (*it).begin();
    iterator_container.push_back(temp_iterator);
    }

    //3- filling the data value
    double data=0.0000;

    //4- Cascade iteration for storing the whole set of combination
    unsigned int K = iterator_container.size();

    while (iterator_container[0] != vector_container[0].end()) {

        //Clear the key_vector before pushing a new key
        std::vector<unsigned int> key_vector;
        key_vector.reserve(variablesTotStates.size());

        //Filling the key_vector 
        for (auto it_key=iterator_container.begin(); it_key!=iterator_container.end(); ++it_key){    
        key_vector.push_back( (*(*it_key)) ); //pushing a new key
        }

        //Assign a new row in the conditional map
        mJointMap.insert( std::make_pair( key_vector, data ) );

    //Increment the odometer by 1 
    //ATTENTION: It must called at the end of the while otherwise the first key get missed...
    ++iterator_container[K-1]; 
    for (int i = K-1; (i > 0) && (iterator_container[i] == vector_container[i].end()); --i) {
        //subtracting the counter
        iterator_container[i] = vector_container[i].begin();
        ++iterator_container[i-1];
    }
    }
    }

    const std::map<std::vector<unsigned int>,double>& JointProbabilityTable::ReturnJointMap(){
    return mJointMap;
    }

    LWSampler::LWSampler(){}

    LWSampler::~LWSampler(){}

    /**
    * It returns a std::pair containing a single sample picked up 
    * from the Bayesian network and the relative weight.
    *
    * @param net the Bayesian network to use for picking up the sample.
    *
    **/
    std::pair< std::vector<unsigned int>, double> LWSampler::ReturnSample(bayonet::Bayesnet& net){

    //Declaring the vector to return
    std::vector<unsigned int> vector_to_return;

    //Declaring the weight
    double weight_to_return = 1.0;

    //Declaring the local variables
    auto topo_list = net.ReturnTopologicalList();
    std::map<unsigned int, std::pair<bool, unsigned int>> sample_map;
    unsigned int tot_nodes = net.ReturnNumberOfNodes();

    //Fill the sample_map and the vector_to_return with zeros
    for(unsigned int i=0; i<tot_nodes; i++){
    auto my_pair = std::make_pair<bool, unsigned int>(false,0);
    auto map_pair = std::make_pair(i, my_pair);
    sample_map.insert(map_pair);
    vector_to_return.push_back(0);
    }

    //Cycle through the topological list
    for(auto it_topo=topo_list.begin(); it_topo!=topo_list.end(); ++it_topo){
    auto in_list = net.ReturnInEdges(*it_topo);
    std::vector<unsigned int> key_vector;
    //Cycle through the in list for creating the key
    for(auto it_in=in_list.begin(); it_in!=in_list.end(); ++it_in){
    auto map_value = sample_map[*it_in]; //return the pair stored inside the map
    bool value_check = map_value.first; //return the boolean of the pair
    unsigned int value_sample = map_value.second; //return the unsigned int of the pair
    if(value_check == false){
        std::cerr << "EXCEPTION: the topological order was not respected!" << std::endl;    
    }else{
        key_vector.push_back(value_sample); //push   
    }
    }
    //Key completed, asking for the sample
    auto sp_to_node = net[*it_topo];
    //Checking if the node is Evidence
    if(sp_to_node.IsEvidence() == true){
    //unsigned int node_evidence = sp_to_node->GetEvidence();
    std::pair<bool,unsigned int> pair_to_store = std::make_pair<bool, unsigned int>(true, sp_to_node.GetEvidence());
    sample_map[*it_topo] = pair_to_store;
    weight_to_return = weight_to_return * sp_to_node.conditionalTable.GetProbability(sp_to_node.GetEvidence(), key_vector);
    }else{
    //If it is not an Evidence then storing the sample in the local index vector
    std::pair<bool,unsigned int> pair_to_store = std::make_pair<bool, unsigned int>(true, sp_to_node.conditionalTable.ReturnSample(key_vector));
    sample_map[*it_topo] = pair_to_store;
    }
    }

    //Fill the vector_to_return with the result stored inside the map
    //This cycle restore the nodes order (not the topological order)
    for(unsigned int i=0; i<tot_nodes; i++){
    auto temp_pair = sample_map[i];
    vector_to_return[i] = temp_pair.second; //Storing the sample into the vector to return
    }

    return std::make_pair(vector_to_return, weight_to_return);
    }

    /**
    * It accumulate samples picking up them from the Bayesian network.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    *
    **/
    std::vector<std::pair<std::vector<unsigned int>, double>> LWSampler::AccumulateSamples(Bayesnet& net, unsigned int cycles){
    std::vector<std::pair<std::vector<unsigned int>, double>> vector_to_return;
    for(unsigned int i=0; i<cycles; i++){
    vector_to_return.push_back(ReturnSample(net));
    }
    return vector_to_return;
    }


    /**
    * It prints the result of the sampling.
    * It is possible to do it for different iterations.
    * The format used is the following :
    * INDEX ..... SAMPLE ..... WEIGHT
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    *
    **/
    void LWSampler::PrintSample(bayonet::Bayesnet& net, unsigned int cycles){
    for(unsigned int i=0; i<cycles; i++){
    auto sample_pair = ReturnSample(net);
    auto sample_vector = sample_pair.first;
    auto weight = sample_pair.second;
    std::cout << i+1 << " ..... ";
    for(auto it=sample_vector.begin(); it!=sample_vector.end(); ++it){
    std::cout << *it << " ";
    }
    std::cout << " ..... " << weight << std::endl;
    }
    }

    /**
    * It creates a Joint Probability table starting from the Bayesian network and sampling for
    * the number of iterations specified.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    * @return it return a Joint Probability Table object
    **/
    JointProbabilityTable LWSampler::ReturnJointProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles){

    //0-Declare the JPT
    JointProbabilityTable joint_table(net.ReturnTotalStates());

    //1-reset JPT
    joint_table.ResetProbabilities();

    
    //auto samples_vector = AccumulateSamples(net, cycles);

    //2-Accumulate samples and Add sample and weight to JPT
    for(unsigned int i=0; i<cycles; i++){
    auto sample_pair = ReturnSample(net);
    joint_table.AddToProbability(sample_pair.first, sample_pair.second);
    }

    //3-Add sample and weight to JPT
    //for(auto it_sample=samples_vector.begin(); it_sample!=samples_vector.end(); ++it_sample){
    //joint_table.AddToProbability(it_sample->first, it_sample->second);
    //}

    //4-Normalize JPT
    joint_table.NormalizeProbabilities();

    //5-Return JPT
    return joint_table;
    }

    /**
    * It creates a Marginal Probability table starting from the Bayesian network and sampling for
    * the number of iterations specified.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    * @return it return a Marginal Probability Table object
    **/
    MarginalProbabilityTable LWSampler::ReturnMarginalProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles){

    MarginalProbabilityTable marginal_table(net.ReturnTotalStates());

    marginal_table.ResetProbabilities();

    //2-Accumulate samples and Add sample and weight to JPT
    for(unsigned int i=0; i<cycles; i++){
    auto sample_pair = ReturnSample(net);
    unsigned int var_counter = 0;
    //Iteration through each element of the sample
    //each element is a variable and the value is the state
    for(auto it_sample=sample_pair.first.begin(); it_sample!=sample_pair.first.end(); ++it_sample){
    marginal_table.AddToProbability(var_counter, *it_sample, sample_pair.second);
    var_counter++;
    } 

    }

    marginal_table.NormalizeProbabilities();

    return marginal_table;
    }

    MarginalProbabilityTable::MarginalProbabilityTable(unsigned int totVariables, unsigned int totStates){

    //Local variables
    std::vector<double> states_vector;
    states_vector.reserve(totStates);
    double probability = (double) 1 / (double) totStates;

    //Filling the states_vector
    for(unsigned int i_stat=0; i_stat<totStates; i_stat++){
    states_vector.push_back(probability);
    }

    //Filling the marginal table
    for(unsigned int i_var=0; i_var<totVariables; i_var++){
    marginalTable.push_back(states_vector);
    //std::cout << " --- " << std::endl;
    }
    }

    /**
    * The constructor of the object.
    *
    * @param variablesTotStatesVector a vector of integers representing the total number of states
    * for each variable.
    **/
    MarginalProbabilityTable::MarginalProbabilityTable(std::vector<unsigned int> variablesTotStatesVector){
    //Cyclyng through the input vector
    for(auto it=variablesTotStatesVector.begin(); it!=variablesTotStatesVector.end(); ++it){
    //Filling the states_vector
    std::vector<double> states_vector;
    states_vector.reserve(*it);
    double probability = (double) 1 / (double) *it;
    for(unsigned int i_stat=0; i_stat<*it; i_stat++){
    states_vector.push_back(probability); 
    }
    marginalTable.push_back(states_vector);
    }
    }

    MarginalProbabilityTable::~MarginalProbabilityTable(){};

    /**
    * It set the probabilities associated with a certain variable
    *
    * @param variableIndex the index of the variable
    * @param stateIndex the index of the state
    * @param probability the value to set
    **/
    bool MarginalProbabilityTable::SetProbability(unsigned int variableIndex, unsigned int stateIndex, double probability){
    if(variableIndex > marginalTable.size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return false; 
    }

    if(variableIndex > marginalTable[variableIndex].size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return false; 
    } 

    marginalTable[variableIndex][stateIndex] = probability;
    return true;
    }

    /**
    * It get the probability associated with a certain variable and state
    *
    * @param variableIndex the index of the variable
    * @param stateIndex the index of the state
    * @param probability
    **/
    bool MarginalProbabilityTable::AddToProbability(unsigned int variableIndex, unsigned int stateIndex, double probability){
    if(variableIndex > marginalTable.size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return false; 
    }

    if(stateIndex > marginalTable[variableIndex].size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return false; 
    } 

    marginalTable[variableIndex][stateIndex] += probability;
    return true;
    }

    /**
    * It get the probability associated with a certain variable and state
    *
    * @param variableIndex the index of the variable
    * @param stateIndex the index of the state
    **/
    double MarginalProbabilityTable::GetProbability(unsigned int variableIndex, unsigned int stateIndex){
    if(variableIndex > marginalTable.size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return false; 
    }

    if(variableIndex > marginalTable[variableIndex].size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return false; 
    } 

    return marginalTable[variableIndex][stateIndex];
    }

    /**
    * It returns the state with the highest probability
    * once specified a variable index. If the variable
    * has many states with the same probability, the first
    * of that state is returned.
    *
    * @param variableIndex
    **/
    unsigned int MarginalProbabilityTable::ReturnMostProbableState(unsigned int variableIndex){
    auto it_max = std::max_element(marginalTable[variableIndex].begin(), marginalTable[variableIndex].end());

    unsigned int counter=0;
    for(auto it_state=marginalTable[variableIndex].begin(); it_state!=marginalTable[variableIndex].end(); ++it_state){
    if(it_state == it_max) return counter;
    counter++;
    }

    return 0;
    }

    /**
    * It set the probabilities associated with a certain variable 
    *
    * @param index the index of the variable
    * @param probabilitiesVector the probabilities of the variable
    **/
    bool MarginalProbabilityTable::SetProbabilities(unsigned int index, std::vector<double> probabilitiesVector){
    if(index > marginalTable.size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return false; 
    }
    marginalTable[index] = probabilitiesVector;
    return true;
    }

    /**
    * It get the probabilities associated with a certain variable
    *
    * @param index the index of the variable
    **/
    std::vector<double> MarginalProbabilityTable::GetProbabilities(unsigned int index){
    if(index > marginalTable.size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    std::vector<double> void_vector;
    void_vector.reserve(1);
    return void_vector; 
    }
    return marginalTable[index];
    }

    /**
    * It reset the probabilities, all values equal to zero
    *
    **/
    void MarginalProbabilityTable::ResetProbabilities(double valueToSet){
    //Iterating through the line
    for(auto it_table=marginalTable.begin(); it_table!=marginalTable.end(); ++it_table){
    //Iterating through the elements
    for(auto it_state=it_table->begin(); it_state!=it_table->end(); ++it_state){
    *it_state = valueToSet;
    }
    }
    }

    /**
    * It normalizes all the probabilities
    *
    **/
    void MarginalProbabilityTable::NormalizeProbabilities(){
    //Iterating through the line
    for(auto it_table=marginalTable.begin(); it_table!=marginalTable.end(); ++it_table){
    //Iterating through the elements
    //for accumulating the divisor
    double accumulator=0.0;
    for(auto it_state=it_table->begin(); it_state!=it_table->end(); ++it_state){
    accumulator += *it_state;
    }

    if(accumulator==0){
    //sometimes the accumulator can be equal to zero
    //In this case a division by 0 is not allowed
    //then all the values are set to zero.
    for(auto it_state=it_table->begin(); it_state!=it_table->end(); ++it_state){
        *it_state = 0;
    }
    }else{
    //normalize the row
    for(auto it_state=it_table->begin(); it_state!=it_table->end(); ++it_state){
        *it_state = *it_state / accumulator;
    }
    }
    }
    }

    /**
    * It prints the probabilities associated with each variable
    *
    **/
    void MarginalProbabilityTable::Print(){
    unsigned int var_counter=0;
    for(auto it_var=marginalTable.begin(); it_var!=marginalTable.end(); ++it_var){
    std::cout << "===== VARIABLE: " << var_counter << " =====" << std::endl;
    unsigned int states_counter=0;
    for(auto it_state=it_var->begin(); it_state!=it_var->end(); ++it_state){
    std::cout << "STATE: " << states_counter << " ..... " << *it_state << std::endl;
    states_counter++;
    } 
    var_counter++;
    }
    }


    /**
    * It prints the probabilities associated with each variable
    *
    **/
    void MarginalProbabilityTable::PrintVariable(unsigned int index){
    if(index > marginalTable.size()){
    std::cerr << "ERROR: Marginal Table out of range index." << std::endl;
    return; 
    }
    std::cout << "===== VARIABLE: " << index << " =====" << std::endl;
    unsigned int states_counter=0;
    for(auto it_state=marginalTable[index].begin(); it_state!=marginalTable[index].end(); ++it_state){
    std::cout << "STATE: " << states_counter << " ..... " << *it_state << std::endl;
    states_counter++;
    } 
    }

    MaximumLikelihoodLearning::MaximumLikelihoodLearning(){}

    MaximumLikelihoodLearning::~MaximumLikelihoodLearning(){}


    /**
    * It uses the training data to update the Conditional Probability
    * Tables of each node in the network.
    *
    * @param net the Bayesian network to update
    * @param trainingDataset a reference to a vector of vector (matrix)
    * containing the state of each node obtained during the training phase.
    * Each line of the matrix must have the same number of values of the 
    * network given as input.
    * @return it return the updated network
    *
    **/
    Bayesnet MaximumLikelihoodLearning::ReturnUpdatedNetwork(Bayesnet net, std::vector<std::vector<unsigned int>>& trainingDataset){

    //1- For each node, Reset to zero the CPTs
    //There are states that are never learnt during
    //the training phase, in this case it is better
    //to initialize the value of the conditionalTable to 1.0
    //and add the dataset sample to them.
    //"Artificial Intelligence: A Modern Approach." chapter 17
    unsigned int tot_nodes = net.ReturnNumberOfNodes();
    for(unsigned int i=0; i<tot_nodes; i++){
    net[i].conditionalTable.ResetProbabilities(1.0); //<--- TRICK from the book 
    }

    //2- For each training vecotr, For each node, take the conditional table, 
    //create a key given the input edges, add +1 to the entry in the table.
    for(auto it_data = trainingDataset.begin(); it_data!=trainingDataset.end(); ++it_data){
    auto sample_vector = *it_data;
    //Iterating all the nodes
    for(unsigned int i=0; i<tot_nodes; i++){
    unsigned int node_state = sample_vector[i];
    auto in_list = net.ReturnInEdges(i);
    std::vector<unsigned int> key_vector;
    key_vector.reserve(in_list.size());
    //Creating the key
    for(auto it_in=in_list.begin(); it_in!=in_list.end(); ++it_in){
        key_vector.push_back(sample_vector[*it_in]);
    }
    //Key vector complete, adding the +1 to the probability 
    net[i].conditionalTable.AddToProbability(node_state, key_vector, 1.0);
    }
    }

    //3- For each node, normalize the CPT
    for(unsigned int i=0; i<tot_nodes; i++){
    net[i].conditionalTable.NormalizeProbabilities();
    }

    //4- return the updated network
    return net;
    }

    RejectionSampler::RejectionSampler(){}

    RejectionSampler::~RejectionSampler(){}

    /**
    * It returns a single sample picking up it from the Bayesian network
    *
    * @param net the Bayesian network to use for picking up the sample.
    *
    **/
    std::vector<unsigned int> RejectionSampler::ReturnSample(bayonet::Bayesnet& net){

    //Declaring the vector to return
    std::vector<unsigned int> vector_to_return;

    //Declaring the local variables
    auto topo_list = net.ReturnTopologicalList();
    std::map<unsigned int, std::pair<bool, unsigned int>> sample_map;
    unsigned int tot_nodes = net.ReturnNumberOfNodes();

    //Fill the sample_map and the vector_to_return with zeros
    for(unsigned int i=0; i<tot_nodes; i++){
    auto my_pair = std::make_pair<bool, unsigned int>(false,0);
    auto map_pair = std::make_pair(i, my_pair);
    sample_map.insert(map_pair);
    vector_to_return.push_back(0);
    }

    //Cycle through the topological list
    for(auto it_topo=topo_list.begin(); it_topo!=topo_list.end(); ++it_topo){
    auto in_list = net.ReturnInEdges(*it_topo);
    std::vector<unsigned int> key_vector;
    //Cycle through the in list for creating the key
    for(auto it_in=in_list.begin(); it_in!=in_list.end(); ++it_in){
    auto map_value = sample_map[*it_in]; //return the pair stored inside the map
    bool value_check = map_value.first; //return the boolean of the pair
    unsigned int value_sample = map_value.second; //return the unsigned int of the pair
    if(value_check == false){
        std::cout << "EXCEPTION: the topological order was not respected!" << std::endl;    
    }else{
        key_vector.push_back(value_sample); //push   
    }
    }
    //Key completed, asking for the sample
    auto sp_to_node = net[*it_topo];
    //Storing the sample in the local index vector
    std::pair<bool,unsigned int> pair_to_store = std::make_pair<bool, unsigned int>(true, sp_to_node.conditionalTable.ReturnSample(key_vector));
    sample_map[*it_topo] = pair_to_store;
    //Storing the sample into the vector to return
    //vector_to_return.push_back(pair_to_store.second);
    }

    //Fill the vector_to_return with the result stored inside the map
    //This cycle restore the nodes order (not the topological order)
    for(unsigned int i=0; i<tot_nodes; i++){
    auto temp_pair = sample_map[i];
    vector_to_return[i] = temp_pair.second; //Storing the sample into the vector to return
    }

    return vector_to_return;
    }

    /**
    * It accumulate samples picking up them from the Bayesian network.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    *
    **/
    std::vector<std::vector<unsigned int>> RejectionSampler::AccumulateSamples(Bayesnet& net, unsigned int cycles){
    std::vector<std::vector<unsigned int>> vector_to_return;
    for(unsigned int i=0; i<cycles; i++){
    vector_to_return.push_back(ReturnSample(net));
    }
    return vector_to_return;
    }

    /**
    * It accumulates samples picking up them from the Bayesian network.
    * It discards the sample that are not coerent with the evidence, the output vector
    * contains only the accepted samples.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    *
    **/
    std::vector<std::vector<unsigned int>> RejectionSampler::AccumulateAndDiscardSamples(Bayesnet& net, unsigned int cycles){
    std::vector<std::vector<unsigned int>> vector_to_return;
    unsigned int tot_nodes = net.ReturnNumberOfNodes();

    //Cycle for accumulating samples
    for(unsigned int i=0; i<cycles; i++){
    auto sample_vector = ReturnSample(net);
    bool discard = false;
    
    //Check all the values and if an evidence is lost it discards the vector
    for(unsigned int nodes_counter=0; nodes_counter<tot_nodes; nodes_counter++){
    if(net[nodes_counter].IsEvidence() == true && net[nodes_counter].GetEvidence() != sample_vector[nodes_counter]){
        discard = true;
        break; //stop the for loop
    }
    }

    //Discard is false, then I push the sample inside the output_vector
    if(discard == false)   vector_to_return.push_back(sample_vector);
    }
    return vector_to_return;
    }

    /**
    * It prints the result of the sampling.
    * It is possible to do it for different iterations.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    *
    **/
    void RejectionSampler::PrintSample(bayonet::Bayesnet& net, unsigned int cycles){

    //Declaring the local variables
    auto topo_list = net.ReturnTopologicalList();
    std::map<unsigned int, std::pair<bool, unsigned int>> sample_map;
    unsigned int tot_nodes = net.ReturnNumberOfNodes();
    unsigned int correct_i_cycle = 0;

    //Cycle for printing
    for (unsigned int i_cycle=0; i_cycle<cycles; i_cycle++){
    correct_i_cycle = i_cycle + 1;
    if(correct_i_cycle <= 9) std::cout << correct_i_cycle << " ...... ";
    else if (correct_i_cycle > 9 && correct_i_cycle <= 99) std::cout << correct_i_cycle << " ..... ";
    else if (correct_i_cycle > 99 && correct_i_cycle <= 999) std::cout << correct_i_cycle << " .... ";
    else std::cout << correct_i_cycle << " ... ";

    //Clear the map at the beginning of the cycle
    sample_map.clear();

    //Fill the outcome vector with zeros
    for(unsigned int i=0; i<tot_nodes; i++){
    auto my_pair = std::make_pair<bool, unsigned int>(false,0);
    auto map_pair = std::make_pair(i, my_pair);
    sample_map.insert(map_pair);
    }

    //Cycle through the topological list
    for(auto it_topo=topo_list.begin(); it_topo!=topo_list.end(); ++it_topo){
    auto in_list = net.ReturnInEdges(*it_topo);
    std::vector<unsigned int> key_vector;
    //Cycle through the in list for creating the key
    for(auto it_in=in_list.begin(); it_in!=in_list.end(); ++it_in){
    auto map_value = sample_map[*it_in]; //return the pair stored inside the map
    bool value_check = map_value.first; //return the boolean of the pair
    unsigned int value_sample = map_value.second; //return the unsigned int of the pair
    if(value_check == false){
        std::cout << "EXCEPTION: the topological order was not respected!" << std::endl;    
    }else{
        key_vector.push_back(value_sample); //push   
    }
    }
    //Key completed, asking for the sample
    auto sp_to_node = net[*it_topo];
    //Storing the sample in the local index vector
    std::pair<bool,unsigned int> pair_to_store = std::make_pair<bool, unsigned int>(true, sp_to_node.conditionalTable.ReturnSample(key_vector));
    sample_map[*it_topo] = pair_to_store;
    //Storing the sample into the vector to return
    //std::cout << pair_to_store.second << " ";
    }

    //Printing the samples in the nodes order (not in the topological one)
    for(auto it_map = sample_map.begin(); it_map!=sample_map.end(); ++it_map){
    auto temp_pair = it_map->second;
    std::cout << temp_pair.second << " ";
    }
    std::cout << std::endl;

    }
    }

    /**
    * It creates a Joint Probability table starting from the Bayesian network and sampling for
    * the number of iterations specified.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    * @return it return a Joint Probability Table object
    **/
    JointProbabilityTable RejectionSampler::ReturnJointProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles){

    //0-Declare the JPT
    JointProbabilityTable joint_table(net.ReturnTotalStates());

    //1-reset JPT
    joint_table.ResetProbabilities();

    //2-Accumulate samples
    auto samples_vector = AccumulateAndDiscardSamples(net, cycles);

    //3-Add sample to JPT
    for(auto it_sample=samples_vector.begin(); it_sample!=samples_vector.end(); ++it_sample){
    joint_table.AddToProbability(*it_sample, 1);
    }

    //4-Normalize JPT
    joint_table.NormalizeProbabilities();

    //5-Return JPT
    return joint_table;
    }

    /**
    * It creates a Marginal Probability table starting from the Bayesian network and sampling for
    * the number of iterations specified.
    *
    * @param net the Bayesian network to use for picking up the sample.
    * @param cycles the number of iterations
    * @return it return a Marginal Probability Table object
    **/
    MarginalProbabilityTable RejectionSampler::ReturnMarginalProbabilityTable(bayonet::Bayesnet& net, unsigned int cycles){

    MarginalProbabilityTable marginal_table(net.ReturnTotalStates());

    marginal_table.ResetProbabilities();

    //2-Accumulate samples and Add sample and weight to JPT
    for(unsigned int i=0; i<cycles; i++){
    auto sample_vector = ReturnSample(net);
    unsigned int var_counter = 0;
    //Iteration through each element of the sample
    //each element is a variable and the value is the state
    for(auto it_sample=sample_vector.begin(); it_sample!=sample_vector.end(); ++it_sample){
    marginal_table.AddToProbability(var_counter, *it_sample, 1.0);
    var_counter++;
    } 

    }

    marginal_table.NormalizeProbabilities();

    return marginal_table;
    }

}


int main()
{
    //Enum is used so that nodes can be named. This makes the next sections easier to read and debug.
    enum node{
    A=0,
    B=1,
    C=2,
    D=3,
    E=4,
    F=5,
    G=6,
    H=7
    };

    std::cout << "This is an implementation of the Variable Elimination method to answer any query for the given Bayesian Network.\n\n";

    //input holds the user query to be calculated
    std::string input = "";

    //Create 3 nodes so max query contains 3 variables. Expanding the code to accomidate more variables is easy but I feel is unnecessary for this project.
    enum node i,j,k;
    std::cout << "Please enter your query by replacing the question mark in Pr(?): ";
    std::cin >> input;
    std::cout << "The answer is:\n";

    //Create a net with 8 nodes where each node can be 1 of 2 states (0 or 1)
    bayonet::Bayesnet myNet({2, 2, 2, 2, 2, 2, 2, 2});

    //It defines the number of iteration to use for the sampler and the sampler to be used.
    unsigned int iterations = 50000;
    bayonet::LWSampler myLWSampler;

    //Create the edges to connect nodes as defined by the figure in the initial problem
    myNet.AddEdge(A, C);
    myNet.AddEdge(B, D);
    myNet.AddEdge(B, E);
    myNet.AddEdge(C, F);
    myNet.AddEdge(D, F);
    myNet.AddEdge(E, H);
    myNet.AddEdge(F, G);
    myNet.AddEdge(F, H);

    //Set the probabilities of each node. True and false are used in place of 0 and 1 for evidences for easier reading.
    myNet[A].conditionalTable.SetProbabilities({}, {0.99000000, 0.01000000});

    myNet[B].conditionalTable.SetProbabilities({}, {0.50000000, 0.50000000});

    myNet[C].conditionalTable.SetProbabilities({true}, {0.95000000, 0.0500});
    myNet[C].conditionalTable.SetProbabilities({false}, {0.99000000, 0.01000000});

    myNet[D].conditionalTable.SetProbabilities({true}, {0.90000000, 0.10000000});
    myNet[D].conditionalTable.SetProbabilities({false}, {0.99000000, 0.01000000});

    myNet[E].conditionalTable.SetProbabilities({true}, {0.40000000, 0.60000000});
    myNet[E].conditionalTable.SetProbabilities({false}, {0.70000000, 0.30000000});

    myNet[F].conditionalTable.SetProbabilities({true, true}, {0.00000000, 1.00000000});
    myNet[F].conditionalTable.SetProbabilities({true, false}, {0.00000000, 1.00000000});
    myNet[F].conditionalTable.SetProbabilities({false, true}, {0.00000000, 1.00000000});
    myNet[F].conditionalTable.SetProbabilities({false, false}, {1.00000000, 0.00000000});

    myNet[G].conditionalTable.SetProbabilities({true}, {0.02000000, 0.98000000});
    myNet[G].conditionalTable.SetProbabilities({false}, {0.95000000, 0.05000000});

    myNet[H].conditionalTable.SetProbabilities({true, true}, {0.10000000, 0.90000000});
    myNet[H].conditionalTable.SetProbabilities({true, false}, {0.20000000, 0.80000000});
    myNet[H].conditionalTable.SetProbabilities({false, true}, {0.30000000, 0.70000000});
    myNet[H].conditionalTable.SetProbabilities({false, false}, {0.90000000, 0.10000000});


    //The output of the sampler is a JointTable, that is a useful way to deal with joint probabilities. 
    bayonet::JointProbabilityTable myJointTable = myLWSampler.ReturnJointProbabilityTable(myNet,iterations);
        
    //If input is 1 variable, set the first node to that variable.
    if (input.length() == 1) //ex. Pr(A)
    {
        switch (input[0]) {
            case 'A': i = A; break;
            case 'B': i = B; break;
            case 'C': i = C; break;
            case 'D': i = D; break;
            case 'E': i = E; break;
            case 'F': i = F; break;
            case 'G': i = G; break;
            case 'H': i = H; break;
            default: std::cout << "ERROR: Capitals Only!\n"; return -1;
        }
        //Print the probability table for the given input.
        std::cout << input << "     Pr(" << input << ")\n";
        myJointTable.PrintMarginal(i, 1, "");
    }
    //If input is 2 variables, set first and second nodes.
    else if (input.length() == 3) //ex. Pr(A|B)
    {
        switch (input[0]) {
            case 'A': i = A; break;
            case 'B': i = B; break;
            case 'C': i = C; break;
            case 'D': i = D; break;
            case 'E': i = E; break;
            case 'F': i = F; break;
            case 'G': i = G; break;
            case 'H': i = H; break;
            default: std::cout << "ERROR: Capitals Only!\n"; return -1;
        }
        switch (input[2]) {
            case 'A': j = A; break;
            case 'B': j = B; break;
            case 'C': j = C; break;
            case 'D': j = D; break;
            case 'E': j = E; break;
            case 'F': j = F; break;
            case 'G': j = G; break;
            case 'H': j = H; break;
            default: std::cout << "ERROR: Capitals Only!\n"; return -1;
        }
        //Print the probability table for the given input.
        std::cout << input[2] << "     " << input[0] << "     Pr(" << input << ")\n";

        myNet[j].SetEvidence(false); //Evidence is 0
        myJointTable = myLWSampler.ReturnJointProbabilityTable(myNet,iterations);
        myJointTable.PrintMarginal(i, 2, "0");

        myNet[j].SetEvidence(true); //Evidence is 1
        myJointTable = myLWSampler.ReturnJointProbabilityTable(myNet,iterations);
        myJointTable.PrintMarginal(i, 2, "1");
    }
    //If input is 3 variables, set all 3 nodes.
    else if (input.length() == 4) //ex. Pr(A|BC)
    {
        switch (input[0]) {
            case 'A': i = A; break;
            case 'B': i = B; break;
            case 'C': i = C; break;
            case 'D': i = D; break;
            case 'E': i = E; break;
            case 'F': i = F; break;
            case 'G': i = G; break;
            case 'H': i = H; break;
            default: std::cout << "ERROR: Capitals Only!\n"; return -1;
        }
        switch (input[2]) {
            case 'A': j = A; break;
            case 'B': j = B; break;
            case 'C': j = C; break;
            case 'D': j = D; break;
            case 'E': j = E; break;
            case 'F': j = F; break;
            case 'G': j = G; break;
            case 'H': j = H; break;
            default: std::cout << "ERROR: Capitals Only!\n"; return -1;
        }
        switch (input[3]) {
            case 'A': k = A; break;
            case 'B': k = B; break;
            case 'C': k = C; break;
            case 'D': k = D; break;
            case 'E': k = E; break;
            case 'F': k = F; break;
            case 'G': k = G; break;
            case 'H': k = H; break;
            default: std::cout << "ERROR: Capitals Only!\n"; return -1;
        }
        //Print the probability table for the given input.
        std::cout << input[2] << "     " << input[3] << "     " << input[0] << "     Pr(" << input << ")\n";

        myNet[j].SetEvidence(false); //Evidence is 00
        myNet[k].SetEvidence(false); 
        myJointTable = myLWSampler.ReturnJointProbabilityTable(myNet,iterations);
        myJointTable.PrintMarginal(i, 3, "00");

        myNet[j].SetEvidence(false); //Evidence is 01
        myNet[k].SetEvidence(true); 
        myJointTable = myLWSampler.ReturnJointProbabilityTable(myNet,iterations);
        myJointTable.PrintMarginal(i, 3, "01");

        myNet[j].SetEvidence(true); //Evidence is 10
        myNet[k].SetEvidence(false); 
        myJointTable = myLWSampler.ReturnJointProbabilityTable(myNet,iterations);
        myJointTable.PrintMarginal(i, 3, "10");

        myNet[j].SetEvidence(true); //Evidence is 11
        myNet[k].SetEvidence(true); 
        myJointTable = myLWSampler.ReturnJointProbabilityTable(myNet,iterations);
        myJointTable.PrintMarginal(i, 3, "11");
    }
    else
    {
        std::cout << "ERROR\n";
    }
    std::cout << std::endl;

    return 0;
}