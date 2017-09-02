#include <iostream>
#include "../btree.h"
#include "catch.hpp"

TEST_CASE("Construct") {
    BTree<int, int> *tree = new BTree<int, int>();
    REQUIRE(tree);

}

TEST_CASE("INSERT") {
    BTree<int, int> *tree = new BTree<int, int>();
    REQUIRE(tree);

    // SECTION("INSERT") {
    // }
}

TEST_CASE("SEARCH") {
    BTree<int, int> *tree = new BTree<int, int>();
    REQUIRE(tree);
}

TEST_CASE("DELETE") {
    BTree<int, int> *tree = new BTree<int, int>();
    REQUIRE(tree);
}
