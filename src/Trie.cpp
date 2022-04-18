#include "Trie.h"

//#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#define TERMINATOR '@'

//using namespace Rcpp;
int Trie::Get_Links_Position(char base)
{
    /*
  Determine the array position of the given base. 
  0, 1, 2, 3, 4 are @ A C G T respectively.
  base: the char to convert to int. Expects either @, A, C, G or T
  return: an int in the inclusive range 0 - 4
  */
    switch (base)
    {
    case 'A':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 3;
    case 'T':
        return 4;
    case TERMINATOR:
        return 0;
    default:
        return -1;
    }
}

trie_node *
Trie::Initialise_Node(char base)
{
    /*
  Initialise a trie node, which is a struct containing the base, the insertion count, 
  and an array to link the attached trie_nodes.
  base: the base rep to create a node of
  return: a pointer to the created node
  */
    trie_node *this_node = new trie_node;
    this_node->base = base;
    this_node->count = 0;
    this_node->end = NULL;

    for (int i = 0; i < 5; i++)
    {
        this_node->links[i] = NULL;
    }

    return this_node;
}

trie_node *
Trie::Initialise_End_Node(char base, int original_seq_index, int current_seq_index)
{
    /* 
  The end node is a special trie node, containing the index of the sequence 
  which ends at this node in the corresponding barcode or hairpin array
  base: the base of this new node, which will always be the terminator character.
  sequence_index: the index in the associated array of sequence which ends at this node.
  */
    trie_node *this_node = Initialise_Node(base);
    end_node *end = new end_node;

    end->original_seq_index = original_seq_index;
    end->current_seq_index = current_seq_index;
    this_node->end = end;
    return this_node;
}

bool Trie::Base_In_Node(trie_node *node, char base)
{
    /*
  Check if the given node contains a link to the given base.
  Checks for the NULL pointer in the position signified by base.
  node: the node to check existance of base in
  base: the base to check for
  return: true if base exists in node->links
          false otherwise.
  */
    if (Get_Links_Position(base) == -1) {
        return false;
    }
    if (node->links[Get_Links_Position(base)] != NULL)
    {
        return true;
    }
    return false;
}

std::vector<trie_node *> Trie::Get_Valid_Links(trie_node *node) {
	if (node == NULL) {
		return {};
	}

	std::vector<trie_node *> out;

	for (int i = 0; i < 5; i++) {
		if (node->links[i] != NULL) {
			out.push_back(node->links[i]);
		}
	}

	return out;
}

trie_node *
Trie::Add_Node(trie_node *node, char base)
{
    /* 
  Adds a trie node to the given node's internal list of nodes,
  returning a pointer to that new node
  Requires the node to not currently contain a link to a node in the array at index base.
  node: a pointer to the node to insert at
  base: the char rep of the node to create
  return: a pointer to the new node created
  */
    node->count++;
    trie_node *new_node = Initialise_Node(base);
    node->links[Get_Links_Position(base)] = new_node;
    return new_node;
}

trie_node *
Trie::Add_End_Node(trie_node *node, char base, int original_seq_index, int current_seq_index)
{
    /* 
  Adds an end node to the trie, which contains all of the 
  same data as a regular node, but with an additional array
  to store the index of the completed hairpin sequences
  node: a pointer to the node to insert at
  base: the char rep of the node to create
  sequence_index: the index of the sequence ending at this node
  return: a pointer to the new node created
  */
    node->count++;
    trie_node *new_node = Initialise_End_Node(base, original_seq_index, current_seq_index);
    node->links[Get_Links_Position(base)] = new_node;
    return new_node;
}

void Trie::Add_String(std::string sequence, int original_seq_index, int current_seq_index) 
{
    trie_node* current_node = head;
    char insert_base;
    current_node = head;


    // loop through each character in the barcode to insert it into the trie
    for (int insert_i = 0; insert_i < (int)sequence.length(); insert_i++)
    {
        insert_base = sequence[insert_i];

        if (Trie::Base_In_Node(current_node, insert_base))
        {
            // if the base is in the current node, simple increment the count
            // and move onto the linked node
            current_node->count++;
            current_node = current_node->links[Trie::Get_Links_Position(insert_base)];
        }
        else
        {
            current_node = Trie::Add_Node(current_node, insert_base);
        }
    }
    // insert the final TERMINATOR @
    // Barcodes aren't always unique, as through using paired reads or dual indexing,
    // the barcode struct will be unique but the barcode sequence or sequenceRev or sequence2 won't necessairly be.
    // so, only add a terminator node if one doesn't always exist.
    if (!Trie::Base_In_Node(current_node, TERMINATOR))
    {
        current_node = Trie::Add_End_Node(current_node, TERMINATOR, original_seq_index, current_seq_index);
    }
    else
    {
        current_node = current_node->links[Trie::Get_Links_Position(TERMINATOR)];
    }
}

// Functions for locating barcodes and hairpins using a Trie initially, and reverting to interative mismatch search if the Trie yields no results
int Trie::Locate_Seq_At_Pos(std::string read, int index, int barcode_length) {
    char base;
    trie_node *current_node;
    end_node *end;
    current_node = head;
    if ((int)read.length() < index) {
        //Rprintf("Short read: %s. Index: %d\n", read.c_str(), index);
        //exit(0);
        return -1;
    }

    if (!current_node) {
        return -1;
    }

    // search from index until we find a TERMINATOR
    for (int j = index; j < index + barcode_length; j++)
    {
        if (j >= (int)read.length()) break;
        base = read[j];
        if (Trie::Base_In_Node(current_node, TERMINATOR))
        {
            // IF the current node can be terminated, return the found hairpin.
            current_node = current_node->links[Trie::Get_Links_Position(TERMINATOR)];
            end = current_node->end;

            return end->original_seq_index;
        }
        else if (Trie::Base_In_Node(current_node, base))
        {
            // If we can continue traversing the trie, move to the next node
            current_node = current_node->links[Trie::Get_Links_Position(base)];
        }
        else
        {
            break;
        }
    }
    // last check if we've reached the end of the read, and the TERMINATOR node exists
    if (Trie::Base_In_Node(current_node, TERMINATOR))
    {
        current_node = current_node->links[Trie::Get_Links_Position(TERMINATOR)];
        end = current_node->end;

        return end->original_seq_index;
    }
    return -1;
}

std::vector<MismatchResult> Trie::Locate_Seq_Mismatches(std::string read, int index, int barcode_length) const {
	std::vector<MismatchResult> out;
	this->SeqMismatchAux(out, read, this->head, index, barcode_length, -1);

	return out;
}

/*
Recursively traverse the trie for exact matches (and mismatches)
against the given string.
*/
void Trie::SeqMismatchAux(std::vector<MismatchResult> &found, const std::string &read, trie_node *current, int index, int left, int mismatchPos) const {
	// if we're at the end of the sequence to check
	if (left == 0) {
		if (this->Base_In_Node(current, TERMINATOR)) {
			MismatchResult n = {current->links[this->Get_Links_Position(TERMINATOR)]->end->original_seq_index, mismatchPos};
			found.push_back(n);
		} 
		return;
	}

	char base = read[index];

	// if base is N, we can match to any barcode character
	// so simply DFS through the trie again 
	// matching a N character does not 'use' up a mismatch
	if (base == 'N') {
		std::vector<trie_node *> valid_nodes = this->Get_Valid_Links(current);
		for (int i = 0; i < valid_nodes.size(); i++) {
			this->SeqMismatchAux(found, read, valid_nodes[i], index + 1, left - 1, mismatchPos);
		}
	} else { // base is not N so check exact matching, then mismatch if we can.
		if (this->Base_In_Node(current, base)) {
			this->SeqMismatchAux(found, read, current->links[this->Get_Links_Position(base)], index + 1, left - 1, mismatchPos);
		} else if (mismatchPos == -1) { // mismatching is only allowed if we can't exact match. All exact matches should be used.
			// mismatchPos == -1 if we haven't yet used up the mismatch
			std::vector<trie_node *> valid_nodes = this->Get_Valid_Links(current);
			for (int i = 0; i < valid_nodes.size(); i++) {
				this->SeqMismatchAux(found, read, valid_nodes[i], index + 1, left - 1, index);
			}
		}
	}
}

int Trie::Locate_Seq_Subsection(std::string read, int section_start, int section_end, int *found_position)
{
    /*
  Search through this read until we locate a known sequence in the given trie. Return that sequences' index.
  Otherwise, if we reach the end of the read, return -1
  trie_head: the head of the trie to search, can be either a barcode trie or a hairpin trie.
  read: the fastq read to search through
  found_position: a pointer to an int which will be used as an out parameter, to store the position in the read the sequence was found at
  return: the original index of the read in the hairpins or barcode array, or -1 if not found
  */
    int i, j;
    char base;
    trie_node *current_node;
    end_node *end;
    for (i = section_start; i < section_end; i++)
    {
        current_node = head;
        // search from i until we find a TERMINATOR
        for (j = i; j < (int)read.length(); j++)
        {
            base = read[j];
            if (Trie::Base_In_Node(current_node, TERMINATOR))
            {
                // IF the current node can be terminated, return the found sequence.
                current_node = current_node->links[Trie::Get_Links_Position(TERMINATOR)];
                end = current_node->end;
                *found_position = i;
                return end->original_seq_index;
            }
            else if (Trie::Base_In_Node(current_node, base))
            {
                // If we can continue traversing the trie, move to the next node
                current_node = current_node->links[Get_Links_Position(base)];
            }
            else
            {
                // else, we have to start searching from the next position in the read.
                break;
            }
        }
        // last check if we've reached the end of the read, and the TERMINATOR node exists
        if (Trie::Base_In_Node(current_node, TERMINATOR))
        {
            current_node = current_node->links[Trie::Get_Links_Position(TERMINATOR)];
            end = current_node->end;
            *found_position = i;

            return end->original_seq_index;
        }
    }
    *found_position = -1;
    return -1;
}

void Trie::Clear_Trie() {
    clear_trie_rec(head);
    head = Initialise_Node('\0');
}
void Trie::clear_trie_rec(trie_node *node)
{
    /* 
  Recursive function to clear a trie. 
  Calls this function on each linked node, then frees this node

  node: a pointer to the current node to free
  */
    int i;
    if (node->end != NULL)
    {
        delete node->end;
    }
    for (i = 0; i < 5; i++)
    {
        if (node->links[i] != NULL)
        {
            clear_trie_rec(node->links[i]);
        }
    }
    delete node;
}
