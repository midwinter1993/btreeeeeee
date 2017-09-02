#ifndef __BNODE_H__
#define __BNODE_H__

#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>


template<typename K, typename V>
class BTree;

template<typename K, typename V>
class BNode {
public:
    typedef K key_t;
    typedef V value_t;
    typedef std::pair<key_t, value_t> kv_t;

    static const size_t ORDER = 2;
    static const size_t MIN_NR_KEY = ORDER - 1;
    static const size_t MAX_NR_KEY = ORDER * 2 - 1;
    static const size_t NPOS = std::numeric_limits<size_t>::max();

    BNode();

private:
    bool isLeaf() const;
    size_t keySize() const;
    bool isFull() const;

    const kv_t& getKV(size_t pos) const;

    const key_t& getKey(size_t pos) const;
    bool keyEqual(size_t pos, const key_t &k) const;

    const value_t& getValue(size_t pos) const;
    void setValue(size_t pos, const value_t &v);

    BNode* getLeftChild(size_t pos) const;
    BNode* getRightChild(size_t pos) const;

private:
    // Search
    std::pair<const BNode*, size_t> search(const key_t &k) const;
    size_t findKey(const key_t &k) const;

    // Insert
    std::pair<kv_t, BNode*> split();
    void insertNotFull(const key_t &k, const value_t &v);

    // Delete
    kv_t findMax() const;
    kv_t findMin() const;
    //
    // Kpos-1     Kpos       Kpos+1
    //        Cpos    Cpos+1
    //
    BNode* mergeChildren(size_t pos);
    void rotate(size_t pos, bool clockwise);

    size_t findSiblingWithMoreKey(size_t pos);
    size_t findSiblingToMerge(size_t pos);

    void remove(const key_t &k);

private:
    friend class BTree<K, V>;

    std::vector<kv_t> kv_pairs_;
    std::vector<BNode*> children_;
    bool is_leaf_;
};

template<typename K, typename V>
BNode<K, V>::BNode(): is_leaf_(true) {
    kv_pairs_.reserve(MAX_NR_KEY);
    children_.reserve(MAX_NR_KEY+1);
}

//
// ---------------------------------------------------------------------------
//
template<typename K, typename V>
bool BNode<K, V>::isLeaf() const {
    return is_leaf_;
}

template<typename K, typename V>
size_t BNode<K, V>::keySize() const {
    return kv_pairs_.size();
}

template<typename K, typename V>
bool BNode<K, V>::isFull() const {
    return keySize() == MAX_NR_KEY;
}

template<typename K, typename V>
const typename BNode<K, V>::kv_t& BNode<K, V>::getKV(size_t pos) const {
    assert(pos < keySize());
    return kv_pairs_[pos];
}

template<typename K, typename V>
const typename BNode<K, V>::key_t& BNode<K, V>::getKey(size_t pos) const {
    return getKV(pos).first;
}

template<typename K, typename V>
bool BNode<K, V>::keyEqual(size_t pos, const key_t &k) const {
    if (keySize() == 0) {
        assert(pos == 0);
        return false;
    } else if (keySize() <= pos) {
        return false;
    }
    return getKey(pos) == k;
}

template<typename K, typename V>
const typename BNode<K, V>::value_t& BNode<K, V>::getValue(size_t pos) const {
    return getKV(pos).second;
}

template<typename K, typename V>
void BNode<K, V>::setValue(size_t pos, const value_t &v) {
    assert(pos < keySize());
    kv_pairs_[pos].second = v;
}

template<typename K, typename V>
BNode<K, V>* BNode<K, V>::getLeftChild(size_t pos) const {
    assert(pos < keySize() + 1 && !isLeaf());
    return children_[pos];
}

template<typename K, typename V>
BNode<K, V>* BNode<K, V>::getRightChild(size_t pos) const {
    assert(pos+1 < keySize() + 1 && !isLeaf());
    return children_[pos+1];
}

//
// ---------------------------------------------------------------------------
// Search
// ---------------------------------------------------------------------------
//
template<typename K, typename V>
size_t BNode<K, V>::findKey(const key_t &k) const {
    auto it = std::lower_bound(kv_pairs_.begin(), kv_pairs_.end(), k,
        [](const kv_t &kv, const key_t &k) -> bool {
            return kv.first < k;
        }
    );
    size_t pos = it - kv_pairs_.begin();
    return pos;
}

template<typename K, typename V>
std::pair<const BNode<K, V>*, size_t> BNode<K, V>::search(const key_t &k) const {
    size_t pos = findKey(k);
    if (keyEqual(pos, k)) {
        return std::make_pair(this, pos);
    }
    if (isLeaf()) {
        return std::make_pair(nullptr, 0);
    }
    auto child = getLeftChild(pos);
    return child->search(k);
}

//
// ---------------------------------------------------------------------------
// Insert
// ---------------------------------------------------------------------------
//
template<typename K, typename V>
std::pair<typename BNode<K, V>::kv_t,
          BNode<K, V>*> BNode<K, V>::split() {
    assert(isFull());
    BNode<K, V> *new_node = new BNode<K, V>();

    new_node->is_leaf_ = is_leaf_;
    std::copy(kv_pairs_.begin() + ORDER,
              kv_pairs_.end(),
              std::back_inserter(new_node->kv_pairs_));
    std::copy(children_.begin() + ORDER,
              children_.end(),
              std::back_inserter(new_node->children_));

    auto ret = std::make_pair(getKV(ORDER-1), new_node);

    kv_pairs_.erase(kv_pairs_.begin() + ORDER-1);
    if (!isLeaf()) {
        children_.erase(children_.begin() + ORDER);
    }
    return ret;
}

template<typename K, typename V>
void BNode<K, V>::insertNotFull(const key_t &k, const value_t &v) {
    assert(!isFull());

    size_t pos = findKey(k);
    if (keyEqual(pos, k)) {
        setValue(pos, v);
    } else if (isLeaf()) {
        kv_pairs_.insert(kv_pairs_.begin() + pos, std::make_pair(k, v));
    } else {
        auto child = getLeftChild(pos);
        if (child->isFull()) {
            auto ret_split = child->split();
            auto &new_kv_pair = ret_split.first;
            auto *new_child = ret_split.second;

            kv_pairs_.insert(kv_pairs_.begin() + pos, new_kv_pair);
            children_.insert(children_.begin() + pos+1, new_child);

            if (k <= new_kv_pair.first) {
                child->insertNotFull(k, v);
            } else {
                new_child->insertNotFull(k, v);
            }
        } else {
            child->insertNotFull(k, v);
        }
    }
}

//
// ---------------------------------------------------------------------------
// Delete
// ---------------------------------------------------------------------------
//
template<typename K, typename V>
typename BNode<K, V>::kv_t BNode<K, V>::findMax() const {
    if (isLeaf()) {
        return kv_pairs_.back();
    }
    return children_.back()->findMax();
}

template<typename K, typename V>
typename BNode<K, V>::kv_t BNode<K, V>::findMin() const {
    if (isLeaf()) {
        return kv_pairs_.front();
    }
    return children_.front()->findMin();
}

template<typename K, typename V>
BNode<K, V>* BNode<K, V>::mergeChildren(size_t pos) {
    //
    // Kpos-1     Kpos       Kpos+1
    //        Cpos    Cpos+1
    //
    auto *left_child = getLeftChild(pos);
    auto *right_child = getRightChild(pos);

    assert(left_child->keySize() <= MIN_NR_KEY &&
           right_child->keySize() <= MIN_NR_KEY);

    left_child->kv_pairs_.push_back(getKV(pos));
    std::copy(right_child->kv_pairs_.begin(),
              right_child->kv_pairs_.end(),
              left_child->kv_pairs_.begin());
    std::copy(right_child->children_.begin(),
              right_child->children_.end(),
              left_child->children_.begin());

    kv_pairs_.erase(kv_pairs_.begin() + pos);
    children_.erase(children_.begin() + pos + 1);

    delete right_child;
    return left_child;
}

template<typename K, typename V>
void BNode<K, V>::rotate(size_t pos, bool clockwise) {
    //
    // Kpos-1     Kpos       Kpos+1
    //        Cpos    Cpos+1
    //
    auto *left_child = getLeftChild(pos);
    auto *right_child = getRightChild(pos);

    if (clockwise) {
        assert(left_child->keySize() > MIN_NR_KEY &&
               right_child->keySize() <= MIN_NR_KEY);

        right_child->kv_pairs_.insert(right_child->kv_pairs_.begin(),
                                      getKV(pos));
        right_child->children_.insert(right_child->children_.begin(),
                                      left_child->children_.back());

        kv_pairs_[pos] = left_child->kv_pairs_.back();

        left_child->kv_pairs_.pop_back();
        left_child->children_.pop_back();
    } else {
        assert(left_child->keySize() <= MIN_NR_KEY &&
               right_child->keySize() > MIN_NR_KEY);

        left_child->kv_pairs_.push_back(getKV(pos));
        left_child->children_.push_back(right_child->children_.front());

        kv_pairs_[pos] = right_child->kv_pairs_.front();

        right_child->kv_pairs_.erase(right_child->kv_pairs_.begin());
        right_child->children_.erase(right_child->children_.begin());
    }
}

template<typename K, typename V>
size_t BNode<K, V>::findSiblingWithMoreKey(size_t pos) {
    if (pos == 0) {
        if (children_[pos+1]->keySize() > MIN_NR_KEY) {
            return pos;
        }
    } else if (pos == children_.size() - 1) {
        if (children_[pos-1]->keySize() > MIN_NR_KEY) {
            return pos-1;
        }
    } else {
        if (children_[pos+1]->keySize() > MIN_NR_KEY) {
            return pos;
        } else if (children_[pos-1]->keySize() > MIN_NR_KEY) {
            return pos-1;
        }
    }
    return NPOS;
}

template<typename K, typename V>
size_t BNode<K, V>::findSiblingToMerge(size_t pos) {
    if (pos == 0) {
        if (children_[pos+1]->keySize() == MIN_NR_KEY) {
            return pos;
        }
    } else if (pos == children_.size() - 1) {
        if (children_[pos-1]->keySize() == MIN_NR_KEY) {
            return pos-1;
        }
    } else {
        if (children_[pos+1]->keySize() == MIN_NR_KEY) {
            return pos;
        } else if (children_[pos-1]->keySize() == MIN_NR_KEY) {
            return pos-1;
        }
    }
    return NPOS;
}

template<typename K, typename V>
void BNode<K, V>::remove(const key_t &k) {
    size_t pos = findKey(k);
    if (isLeaf()) {
        if (keyEqual(pos, k)) {
            kv_pairs_.erase(kv_pairs_.begin() + pos);
        }
    } else {
        if (keyEqual(pos, k)) {
            auto *left_child = getLeftChild(pos);
            auto *right_child = getRightChild(pos);

            if (left_child->keySize() > MIN_NR_KEY) {
                kv_t kv = left_child->findMax();
                kv_pairs_[pos] = kv;
                left_child->remove(kv.first);
            } else if (right_child->keySize() > MIN_NR_KEY) {
                kv_t kv = right_child->findMin();
                kv_pairs_[pos] = kv;
                right_child->remove(kv.first);
            } else {
                mergeChildren(pos);
                left_child->remove(k);
            }
        } else {
            size_t pivot = findSiblingWithMoreKey(pos);
            if (pivot != NPOS) {
                bool clockwise = pivot+1 == pos ? true : false;
                rotate(pivot, clockwise);
                getLeftChild(pos)->remove(k);
            } else {
                size_t parent = findSiblingToMerge(pos);
                assert(parent != NPOS);
                auto *new_node = mergeChildren(parent);
                new_node->remove(k);
            }
        }
    }
}


#endif /* ifndef __BNODE_H__ */
