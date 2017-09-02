#ifndef __BTREE_H__
#define __BTREE_H__

#include "bnode.h"

#include <cassert>
#include <iostream>


template<typename K, typename V>
class BTree {
public:
    typedef K key_t;
    typedef V value_t;

    BTree();
    void insert(const key_t &k, const value_t &v);
    value_t search(const key_t &k) const;
    void remove(const key_t &k);
private:
    BNode<K, V>* makeRootNode(bool is_leaf) const;
    // BNode<K, V>* makeLeafNode() const;
    // BNode<K, V>* makeInternalNode() const;

private:
    BNode<K, V> *root_;
};

template<typename K, typename V>
BTree<K, V>::BTree(): root_(makeRootNode(true)) {
}

template<typename K, typename V>
BNode<K, V>* BTree<K, V>::makeRootNode(bool is_leaf) const {
    auto *root = new BNode<K, V>();
    root->is_leaf_ = is_leaf;
    return root;
}

template<typename K, typename V>
void BTree<K, V>::insert(const key_t &k, const value_t &v) {
    if (root_->isFull()) {
        auto *new_root = makeRootNode(false);
        new_root->children_.push_back(root_);
        root_ = new_root;

        root_->insertNotFull(k, v);
    } else {
        root_->insertNotFull(k, v);
    }
}

template<typename K, typename V>
typename BTree<K, V>::value_t BTree<K, V>::search(const key_t &k) const {
    auto ret = root_->search(k);
    auto *node = ret.first;
    auto pos = ret.second;
    if (!node) {
        return value_t{};
    } else {
        return node->getValue(pos);
    }
}

template<typename K, typename V>
void BTree<K, V>::remove(const key_t &k) {
    root_->remove(k);
}

#endif /* ifndef __BTREE_H__ */
