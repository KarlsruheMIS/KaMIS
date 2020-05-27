/******************************************************************************
 * bucket_pq.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef BUCKET_PQ_EM8YJPA9
#define BUCKET_PQ_EM8YJPA9

#include <limits>
#include <unordered_map>
#include <vector>
#include <assert.h>


class bucket_pq {
        public:
                bucket_pq( const int max_prio ); 

                virtual ~bucket_pq() {};

                int size();  
                void insert(int id, int gain); 
                bool empty();

                int maxValue();
                int maxElement();
                int deleteMax();

                void decreaseKey(int node, int newGain);
                void increaseKey(int node, int newGain);

                void changeKey(int element, int newKey);
                int getKey(int element);
                void deleteNode(int node);

                bool contains(int node);
        private:
                int     m_elements;
                unsigned   m_max_idx; //points to the non-empty bucket with the largest priority
                
                std::unordered_map<int, std::pair<int, int> > m_queue_index;
                std::vector<std::vector<int>> m_buckets;
};

inline bucket_pq::bucket_pq( const int max_prio ) {
        m_elements  = 0;
        m_max_idx   = 0;

        m_buckets.resize(max_prio + 1);
}

inline int bucket_pq::size() {
        return m_elements;  
}

inline void bucket_pq::insert(int node, int prio) {
        unsigned address = prio;
        if(address > m_max_idx) {
                m_max_idx = address; 
        }
       
        m_buckets[address].push_back( node ); 
        m_queue_index[node].first  = m_buckets[address].size() - 1; //store position
        m_queue_index[node].second = prio;

        m_elements++;
}

inline bool bucket_pq::empty( ) {
        return m_elements == 0;        
}

inline int bucket_pq::maxValue( ) {
        return m_max_idx;        
}

inline int bucket_pq::maxElement( ) {
        return m_buckets[m_max_idx].back();        
}

inline int bucket_pq::deleteMax() {
       int node = m_buckets[m_max_idx].back();
       m_buckets[m_max_idx].pop_back();
       m_queue_index.erase(node);

       if( m_buckets[m_max_idx].size() == 0 ) {
             //update max_idx
             while( m_max_idx != 0 )  {
                     m_max_idx--;
                     if(m_buckets[m_max_idx].size() > 0) {
                        break;
                     }
             }
       }

       m_elements--;
       return node;        
}

inline void bucket_pq::decreaseKey(int node, int new_prio) {
        changeKey( node, new_prio );
}

inline void bucket_pq::increaseKey(int node, int new_prio) {
        changeKey( node, new_prio );
}

inline int bucket_pq::getKey(int node) {
        return m_queue_index[node].second;        
}

inline void bucket_pq::changeKey(int node, int new_prio) {
        deleteNode(node);
        insert(node, new_prio);
}

inline void bucket_pq::deleteNode(int node) {
        assert(m_queue_index.find(node) != m_queue_index.end());
        int in_bucket_idx = m_queue_index[node].first;
        int  old_gain      = m_queue_index[node].second;
        unsigned address    = old_gain;

        if( m_buckets[address].size() > 1 ) {
                //swap current element with last element and pop_back
                m_queue_index[m_buckets[address].back()].first = in_bucket_idx; // update helper structure
                std::swap(m_buckets[address][in_bucket_idx], m_buckets[address].back());
                m_buckets[address].pop_back();
        } else {
                //size is 1
                m_buckets[address].pop_back();
                if( address == m_max_idx ) {
                        //update max_idx
                        while( m_max_idx != 0 )  {
                                m_max_idx--;
                                if(m_buckets[m_max_idx].size() > 0) {
                                        break;
                                }
                        }

                }
        }

        m_elements--;
        m_queue_index.erase(node);
}

inline bool bucket_pq::contains(int node) {
        return m_queue_index.find(node) != m_queue_index.end(); 
}


#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */
