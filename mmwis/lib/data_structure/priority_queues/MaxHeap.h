//
// Created by alex on 09.02.20.
//

#ifndef COMPONENTS_MAXHEAP_H
#define COMPONENTS_MAXHEAP_H

#include <limits>
#include <vector>
#include <unordered_map>
#include <execinfo.h>
#include <definitions.h>

template < typename Data, typename Key >
class QElem {
public:
    QElem( Data data, Key key, int index ) : m_data(data), m_key (key), m_index(index) {};
    virtual ~QElem() {};

    Data & get_data() {
        return m_data;
    }

    void set_data(Data & data) {
        m_data = data;
    }

    Key get_key() {
        return m_key;
    }

    void set_key(Key key) {
        m_key = key;
    }

    int get_index() {
        return m_index;
    }

    void set_index(int index) {
        m_index = index;
    }

private:
    Data m_data;
    Key  m_key;
    int  m_index; // the index of the element in the heap

};

template<typename Key>
class MaxHeap {
public:

    struct Data {
        NodeID node;
        Data( NodeID node ) : node(node) {};
    };

    typedef QElem<Data, Key> PQElement;

    MaxHeap() {};
    ~MaxHeap() = default;

    NodeID size();
    bool empty();

    bool contains(NodeID node);
    void insert(NodeID id, Key key);

    NodeID deleteMax();
    void deleteNode(NodeID node);
    NodeID maxElement();
    Key maxValue();

    void decreaseKey(NodeID node, Key key);
    void increaseKey(NodeID node, Key key);
    void changeKey(NodeID node, Key key);
    Key getKey(NodeID node);

    void clear();

private:
    std::vector< PQElement >               m_elements; // elements that contain the data
    std::unordered_map<NodeID, int>   m_element_index; // stores index of the node in the m_elements array
    std::vector< std::pair<Key, int> >     m_heap;     // key and index in elements (pointer)

    void siftUp( int pos );
    void siftDown( int pos );

};

template<typename Key>
inline Key MaxHeap<Key>::maxValue() {
    return m_heap[0].first;
};

template<typename Key>
inline NodeID MaxHeap<Key>::maxElement() {
    return m_elements[m_heap[0].second].get_data().node;
};

template<typename Key>
inline void MaxHeap<Key>::siftDown( int pos ) {

    int curKey   = m_heap[pos].first;
    int lhsChild = 2*pos+1;
    int rhsChild = 2*pos+2;
    if( rhsChild < (int) m_heap.size() ) {

        int lhsKey = m_heap[lhsChild].first;
        int rhsKey = m_heap[rhsChild].first;

        if( lhsKey < curKey && rhsKey < curKey) {
            return; // we are done
        } else {
            //exchange with the larger one (maxHeap)
            int swap_pos = lhsKey > rhsKey ? lhsChild : rhsChild;
            std::swap( m_heap[pos], m_heap[swap_pos]);

            int element_pos = m_heap[pos].second;
            m_elements[element_pos].set_index(pos);

            element_pos = m_heap[swap_pos].second;
            m_elements[element_pos].set_index(swap_pos);

            siftDown(swap_pos);
            return;
        }

    } else if ( lhsChild < (int)m_heap.size()) {
        if( m_heap[pos].first < m_heap[lhsChild].first) {
            std::swap( m_heap[pos], m_heap[lhsChild]);

            int element_pos = m_heap[pos].second;
            m_elements[element_pos].set_index(pos);

            element_pos = m_heap[lhsChild].second;
            m_elements[element_pos].set_index(lhsChild);

            siftDown(lhsChild);
            return;
        } else {
            return; // we are done
        }
    }
}

template<typename Key>
inline void MaxHeap<Key>::siftUp( int pos ) {
    if( pos > 0 ) {
        int parentPos = (int)(pos-1)/2;
        if(  m_heap[parentPos].first < m_heap[pos].first) {
            //heap condition not fulfulled
            std::swap(m_heap[parentPos], m_heap[pos]);

            int element_pos = m_heap[pos].second;
            m_elements[element_pos].set_index(pos);

            // update the heap index in the element
            element_pos = m_heap[parentPos].second;
            m_elements[element_pos].set_index(parentPos);

            siftUp( parentPos );
        }

    }
}

template<typename Key>
inline NodeID MaxHeap<Key>::size() {
    return m_heap.size();
}

template<typename Key>
inline bool MaxHeap<Key>::empty( ) {
    return m_heap.empty();
}

template<typename Key>
inline void MaxHeap<Key>::insert(NodeID node, Key key) {
    if( m_element_index.find(node) == m_element_index.end() ) {
        int element_index =  m_elements.size();
        int heap_size     =  m_heap.size();

        m_elements.push_back( PQElement(Data(node), key, heap_size) );
        m_heap.push_back( std::pair< Key, int>(key, element_index) );
        m_element_index[node] = element_index;
        siftUp( heap_size );
    }
}

template<typename Key>
inline void MaxHeap<Key>::deleteNode(NodeID node) {
    int element_index = m_element_index[node];
    int heap_index    = m_elements[element_index].get_index();

    m_element_index.erase(node);

    std::swap( m_heap[heap_index], m_heap[m_heap.size() - 1]);
    //update the position of its element in the element array
    m_elements[m_heap[heap_index].second].set_index(heap_index);

    // we dont want holes in the elements array -- delete the deleted element from the array
    if(element_index != (int)(m_elements.size() - 1)) {
        std::swap( m_elements[element_index], m_elements[m_elements.size() - 1]);
        m_heap[ m_elements[element_index].get_index() ].second = element_index;
        int cnode              = m_elements[element_index].get_data().node;
        m_element_index[cnode] = element_index;
    }

    m_elements.pop_back();
    m_heap.pop_back();

    if( m_heap.size() > 1 && heap_index < (int)m_heap.size() ) {
        //fix the max heap property
        siftDown(heap_index);
        siftUp(heap_index);
    }
};

template<typename Key>
inline NodeID MaxHeap<Key>::deleteMax() {
    if( m_heap.size() > 0) {
        int element_index = m_heap[0].second;
        int node = m_elements[element_index].get_data().node;
        m_element_index.erase(node);

        m_heap[0] = m_heap[m_heap.size() - 1];
        //update the position of its element in the element array
        m_elements[m_heap[0].second].set_index(0);

        // we dont want holes in the elements array -- delete the deleted element from the array
        if(element_index != (int)(m_elements.size() - 1)) {
            m_elements[element_index] = m_elements[m_elements.size() - 1];
            m_heap[ m_elements[element_index].get_index() ].second = element_index;
            int cnode              = m_elements[element_index].get_data().node;
            m_element_index[cnode] = element_index;
        }

        m_elements.pop_back();
        m_heap.pop_back();

        if( m_heap.size() > 1) {
            //fix the heap property
            siftDown(0);
        }

        return node;
    }

    return -1;
}

template<typename Key>
inline void MaxHeap<Key>::changeKey(NodeID node, Key key) {
    Key old_key = m_heap[m_elements[m_element_index[node]].get_index()].first;
    if(old_key > key ) {
        decreaseKey(node, key);
    } else if (old_key < key ) {
        increaseKey(node, key);
    }
};

template<typename Key>
inline void MaxHeap<Key>::decreaseKey(NodeID node, Key key) {
    ASSERT_TRUE(m_element_index.find(node) != m_element_index.end());
    int queue_idx = m_element_index[node];
    int heap_idx  = m_elements[queue_idx].get_index();
    m_elements[queue_idx].set_key(key);
    m_heap[heap_idx].first = key;
    siftDown(heap_idx);
}

template<typename Key>
inline void MaxHeap<Key>::increaseKey(NodeID node, Key key) {
    ASSERT_TRUE(m_element_index.find(node) != m_element_index.end());
    int queue_idx = m_element_index[node];
    int heap_idx  = m_elements[queue_idx].get_index();
    m_elements[queue_idx].set_key(key);
    m_heap[heap_idx].first = key;
    siftUp(heap_idx);
}

template<typename Key>
inline Key MaxHeap<Key>::getKey(NodeID node) {
    return m_heap[m_elements[m_element_index[node]].get_index()].first;
};


template<typename Key>
inline bool MaxHeap<Key>::contains(NodeID node) {
    return m_element_index.find(node) != m_element_index.end();
}

template<typename Key>
inline void MaxHeap<Key>::clear() {
    m_elements.clear();
    m_element_index.clear();
    m_heap.clear();
}

#endif //COMPONENTS_MAXHEAP_H
