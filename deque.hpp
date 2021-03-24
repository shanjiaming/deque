#ifndef SJTU_DEQUE_HPP
#define SJTU_DEQUE_HPP

#include "exceptions.hpp"

#include <cstddef>

namespace sjtu {

    template<typename T>
    class vector {
    private:
        T **a;
//T**a;//FIXME MemLeak
        size_t space, num;
        static const size_t InitSize = 5;

        void doubleSpace() {
            T **tmp;
//T**tmp;//FIXME MemLeak
            space *= 2;
            tmp = new T *[space];
            for (int i = 0; i < num; ++i)
                tmp[i] = a[i];
            delete[] a;
            a = tmp;
        }

    public:
        /**
         * TODO
         * a type for actions of the elements of a vector, and you should write
         *   a class named const_iterator with same interfaces.
         */
        /**
         * you can see RandomAccessIterator at CppReference for help.
         */
//        class const_iterator;
        class iterator;

        typedef iterator const_iterator;

        class iterator {
            friend vector<T>;
        private:
            const vector<T> *vec;
            int pos;
        public:
            iterator(const vector<T> *v, const int p) : vec(v), pos(p) {}

            /**
             * return a new iterator which pointer n-next elements
             * as well as operator-
             */
            iterator operator+(const int &n) const {
                return iterator(vec, pos + n);
            }

            iterator operator-(const int &n) const {
                return iterator(vec, pos - n);
            }

            // return the distance between two iterators,
            // if these two iterators point to different vectors, throw invaild_iterator.
            int operator-(const iterator &rhs) const {
                if (vec != rhs.vec)
                    throw invalid_iterator();
                return pos - rhs.pos;
            }

            iterator &operator+=(const int &n) {
                pos = pos + n;
                return *this;
            }

            iterator &operator-=(const int &n) {
                pos = pos - n;
                return *this;
            }

            /**
             * TODO iter++
             */
            iterator operator++(int) {
                auto tmp = *this;
                *this += 1;
                return tmp;
            }

            /**
             * TODO ++iter
             */
            iterator &operator++() {
                return *this += 1;
            }

            /**
             * TODO iter--
             */
            iterator operator--(int) {
                auto tmp = *this;
                *this -= 1;
                return tmp;
            }

            /**
             * TODO --iter
             */
            iterator &operator--() {
                return *this -= 1;
            }

            /**
             * TODO *it
             */
            T &operator*() const {
                return *(vec->a[pos]);
            }

            /**
             * a operator to check whether two iterators are same (pointing to the same memory address).
             */
            bool operator==(const iterator &rhs) const {
                return rhs.vec == vec && rhs.pos == pos;
            }
//            bool operator==(const const_iterator &rhs) const {
//                return rhs.vec == vec && rhs.pos == pos;
//            }
            /**
             * some other operator for iterator.
             */
            bool operator!=(const iterator &rhs) const {
                return !(rhs == *this);
            }
//            bool operator!=(const const_iterator &rhs) const {
//                return !(rhs == *this);
//            }
        };

        /**
         * TODO Constructs
         * Atleast two: default constructor, copy constructor
         */
        vector() : space(InitSize), num(0) {
            a = new T *[InitSize];
        }

        vector(const vector &other) : space(other.space), num(other.num) {
            a = new T *[space];
            for (int i = 0; i < num; ++i) {
                a[i] = new T(*other.a[i]);
            }
        }

        /**
         * TODO Destructor
         */
        ~vector() {
            clear();
        }

        /**
         * TODO Assignment operator
         */
        vector &operator=(const vector &other) {
            if (this == &other)
                return *this;
            clear();
            a = new T *[space];
            for (int i = 0; i < num; ++i) {
                a[i] = new T(*other.a[i]);
            }
            return *this;
        }

        /**
         * assigns specified element with bounds checking
         * throw index_out_of_bound if pos is not in [0, space)
         */
        T &at(const size_t &pos) {
            if (pos < 0 || pos >= num)
                throw index_out_of_bound();
            return *a[pos];
        }

        const T &at(const size_t &pos) const {
            if (pos < 0 || pos >= num)
                throw index_out_of_bound();
            return *a[pos];
        }

        /**
         * assigns specified element with bounds checking
         * throw index_out_of_bound if pos is not in [0, space)
         * !!! Pay attentions
         *   In STL this operator does not check the boundary but I want you to do.
         */
        T &operator[](const size_t &pos) {
            if (pos < 0 || pos >= num)
                throw index_out_of_bound();
            return *a[pos];
        }

        const T &operator[](const size_t &pos) const {
            if (pos < 0 || pos >= num)
                throw index_out_of_bound();
            return *a[pos];
        }

        /**
         * access the first element.
         * throw container_is_empty if space == 0
         */
        const T &front() const {
            if (!num)
                throw container_is_empty();
            return *a[0];
        }

        /**
         * access the last element.
         * throw container_is_empty if space == 0
         */
        const T &back() const {
            if (!num)
                throw container_is_empty();
            return *a[num - 1];
        }

        /**
         * returns an iterator to the beginning.
         */
        iterator begin() {
            return iterator(this, 0);
        }

        const_iterator cbegin() const {
            return const_iterator(this, 0);
        }

        /**
         * returns an iterator to the end.
         */
        iterator end() {
            return iterator(this, num);
        }

        const_iterator cend() const {
            return const_iterator(this, num);
        }

        /**
         * checks whether the container is empty
         */
        bool empty() const {
            return !num;
        }

        /**
         * returns the number of elements
         */
        size_t size() const {
            return num;
        }

        /**
         * clears the contents
         */
        void clear() {
            for (int i = 0; i < num; i++)
                delete a[i];
            delete[] a;
            a = nullptr;
        }

        /**
         * inserts value before pos
         * returns an iterator pointing to the inserted value.
         */
        iterator insert(iterator pos, const T &value) {
//            return insert(pos.pos, value);
        }

        /**
         * inserts value at index ind.
         * after inserting, this->at(ind) == value
         * returns an iterator pointing to the inserted value.
         * throw index_out_of_bound if ind > space (in this situation ind can be space because after inserting the space will increase 1.)
         */
        iterator insert(const size_t &ind, const T &value) {
            if (ind < 0 || ind > num)
                throw index_out_of_bound();
            if (num == space)
                doubleSpace();
            for (size_t i = num; i > ind; --i)
                a[i] = a[i - 1];
            a[ind] = new T(value);
            ++num;
            return iterator(this, ind);
        }

        /**
         * removes the element at pos.
         * return an iterator pointing to the following element.
         * If the iterator pos refers the last element, the end() iterator is returned.
         */
        iterator erase(iterator pos) {
//            return erase(pos.pos);
        }

        /**
         * removes the element with index ind.
         * return an iterator pointing to the following element.
         * throw index_out_of_bound if ind >= space
         */
        iterator erase(const size_t &ind) {
            if (ind < 0 || ind >= num)
                throw index_out_of_bound();
            delete a[ind];
            for (size_t i = ind; i < num; ++i)
                a[i] = a[i + 1];
            --num;
            return iterator(this, ind);
        }

        /**
         * adds an element to the end.
         */
        void push_back(const T &value) {
            if (num == space)
                doubleSpace();
            a[num++] = new T(value);
        }

        /**
         * remove the last element from the end.
         * throw container_is_empty if space() == 0
         */
        void pop_back() {
            if (!num)
                throw container_is_empty();
            --num;
            delete a[num];
        }
    };


    template<class T>
    class deque {
    private:

        using Node = T *;

        static const int Nmax = 400;

        typedef int Address;

        struct Block {
            friend deque;
            Address next = -1;      //next means the address of next Block in the file.
            int num = 0;            //num means how many valid nodes are there in this Block.
//        Node nodes[Nmax];       //This means 0 <= num <= Nmax
            vector<Node> nodes;
        };

        vector<Block> blocks;

        Address last_block_address;

        size_t num;

        void getblock(Address x, Block &b) {
            return blocks[x];
        }

        void putblock(Address x, const Block &b) {
            blocks[x] = b;
        }

        void putblockend(const Block &b) {
            blocks.push_back(b);
        }

        int lower_bound(vector<Node> v, Node erase_node) {//FIXME 根本就没有lower_bound逻辑，气死我了。
            int sz = v.size();
            for (int i = 0; i < sz; ++i) {
                if (v[i] == erase_node) {
                    return i;
                }
            }
            return sz;
        }



    public:
        class iterator;

        typedef iterator const_iterator;

//	class const_iterator;
        class iterator {
            friend deque<T>;
        private:
            deque *deque_ptr;
            int block_th;
            int node_th;
            int redundancy;

            iterator(deque *_deque_ptr, int _block_th, int _node_th, int _redundancy = 0) : deque_ptr(_deque_ptr),
                                                                                            block_th(_block_th),
                                                                                            node_th(_node_th),
                                                                                            redundancy(_redundancy) {
            };
        public:
            /**
             * return a new iterator which pointer n-next elements
             *   if there are not enough elements, iterator becomes invalid
             * as well as operator-
             */
            iterator operator+(const int &n) const {
                //TODO
                iterator ret = *this;
                ret += n;
                return ret;
            }

            iterator operator-(const int &n) const {
                //TODO
                iterator ret = *this;
                ret -= n;
                return ret;
            }

            // return th distance between two iterator,
            // if these two iterators points to different vectors, throw invaild_iterator.
            int operator-(const iterator &rhs) const {
                //TODO

            }

            iterator &operator+=(const int &n) {
                //TODO
                if (n < 0)
                    return (*this -= (-n));
                int nchange = n + node_th + redundancy;
                redundancy = 0;
                while (block_th < deque_ptr->blocks.size()) {
                    if (nchange < deque_ptr->blocks[block_th].num) {
                        node_th = nchange;
                        return *this;
                    }
                    nchange -= deque_ptr->blocks[block_th].num;
                    ++block_th;
                }
                redundancy = nchange;
                node_th = 0;
                return *this;
            }

            iterator &operator-=(const int &n) {
                if (n < 0)
                    return (*this += (-n));
                int nchange = -n + node_th + redundancy;
                redundancy = 0;
                while (block_th >= 0) {
                    if (nchange >= 0) {
                        node_th = nchange;
                        return *this;
                    }
                    --block_th;
                    nchange += deque_ptr->blocks[block_th].num;
                }
                redundancy = nchange;
                node_th = 0;
                return *this;
            }

            /**
             * TODO iter++
             */
            iterator operator++(int) {
                iterator ret = *this;
                *this += 1;
                return ret;
            }

            /**
             * TODO ++iter
             */
            iterator &operator++() {
                return *this += 1;
            }

            /**
             * TODO iter--
             */
            iterator operator--(int) {
                iterator ret = *this;
                *this -= 1;
                return ret;
            }

            /**
             * TODO --iter
             */
            iterator &operator--() {
                return *this -= 1;
            }

            /**
             * TODO *it
             * 		throw if iterator is invalid
             */
            T &operator*() const {
                if(redundancy != 0) throw invalid_iterator();
                return deque_ptr->blocks[block_th][node_th];
            }

            /**
             * TODO it->field
             * 		throw if iterator is invalid
             */
            T *operator->() const noexcept {
                return &**this;
            }

            /**
             * a operator to check whether two iterators are same (pointing to the same memory).
             */
            bool operator==(const iterator &rhs) const {}
//		bool operator==(const const_iterator &rhs) const {}
            /**
             * some other operator for iterator.
             */
            bool operator!=(const iterator &rhs) const {}
//		bool operator!=(const const_iterator &rhs) const {}
        };
/*	class const_iterator {
		// it should has similar member method as iterator.
		//  and it should be able to construct from an iterator.
		private:
			// data members.
		public:
			const_iterator() {
				// TODO
			}
			const_iterator(const const_iterator &other) {
				// TODO
			}
			const_iterator(const iterator &other) {
				// TODO
			}
			// And other methods in iterator.
			// And other methods in iterator.
			// And other methods in iterator.
	};*/
        /**
         * TODO Constructors
         */
        deque() {
            last_block_address = -1;
            num = 0;
        }

        deque(const deque &other) {
            blocks = other.blocks;
            last_block_address = other.last_block_address;
            num = other.num;
        }

        /**
         * TODO Deconstructor
         */
        ~deque() {}

        /**
         * TODO assignment operator
         */
        deque &operator=(const deque &other) {
            blocks = other.blocks;
            last_block_address = other.last_block_address;
            num = other.num;
            return *this;
        }

        /**
         * access specified element with bounds checking
         * throw index_out_of_bound if out of bound.
         */
        T &at(const size_t &pos) {
            try{
                auto ret = *(begin() + pos);
                return ret;
            } catch(invalid_iterator) {
                throw index_out_of_bound();
            }
        }

        const T &at(const size_t &pos) const {
            try{
                auto ret = *(cbegin() + pos);
                return ret;
            } catch(invalid_iterator) {
                throw index_out_of_bound();
            }
            //return *(cbegin() + pos);//FIXME 这个版本暂时不能这么写，因为const_iterator是假const。
        }

        T &operator[](const size_t &pos) {
            return at(pos);
        }

        const T &operator[](const size_t &pos) const {
            return at(pos);
        }

        /**
         * access the first element
         * throw container_is_empty when the container is empty.
         */
        const T &front() const {
            if (empty()) throw container_is_empty();
            return *cbegin();
        }

        /**
         * access the last element
         * throw container_is_empty when the container is empty.
         */
        const T &back() const {
            if (empty()) throw container_is_empty();
            return *(--cend());
        }

        /**
         * returns an iterator to the beginning.
         */
        iterator begin() {
            return iterator(this, 0, 0, 0);
        }

        const_iterator cbegin() const {}

        /**
         * returns an iterator to the end.
         */
        iterator end() {
            return iterator(this, last_block_address, blocks[last_block_address].size() - 1, 1);
        }

        const_iterator cend() const {}

        /**
         * checks whether the container is empty.
         */
        bool empty() const {
            return num == 0;
        }

        /**
         * returns the number of elements
         */
        size_t size() const {
            return num;
        }

        /**
         * clears the contents
         */
        void clear() {
            blocks.clear();
            num = 0;
            last_block_address = -1;
        }

        /**
         * inserts elements at the specified locat on in the container.
         * inserts value before pos
         * returns an iterator pointing to the inserted value
         *     throw if the iterator is invalid or it point to a wrong place.
         */
        iterator insert(iterator pos, const T &value) {
            T* poi = new T (value);
            Block bl;
            getblock(0, bl);
            Address block_pos = 0;
            while (bl.num && bl.nodes[bl.num - 1] < insert_node) {
                if (bl.next == -1) break;
                block_pos = bl.next;
                getblock(bl.next, bl);
            }
            int insert_pos = lower_bound(bl.nodes, insert_node);
            if (bl.num < Nmax) {//insert directly
                if (bl.num == insert_pos) {
                    bl.nodes.push_back(insert_node);
//                bl.nodes[insert_pos] = insert_node;
                } else {
                    Node temp = insert_node;
                    Node record;
                    for (int i = insert_pos; i <= bl.num; ++i) {
                        record = bl.nodes[i];
                        bl.nodes[i] = temp;
                        temp = record;
                    }
                }
                ++bl.num;
                putblock(block_pos, bl);
            } else {// insert with split in insert pos
                Block bl_end;
                if (insert_pos != bl.num) {
                    for (int i = insert_pos; i < bl.num; ++i) {
                        bl_end.nodes[i - insert_pos] = bl.nodes[i];
                    }
                    bl.nodes[insert_pos] = insert_node;
                    bl_end.num = bl.num - insert_pos;
                    bl.num = insert_pos + 1;
                } else {

//                bl_end.nodes[0] = insert_node;
                    bl_end.nodes.push_back(insert_node);
                    bl_end.num = bl.num - insert_pos + (insert_pos == bl.num);
                }
                bl_end.next = bl.next;
                last_block_address = bl.next = blocks.size();
                //this is change mark
                putblock(block_pos, bl);
                putblockend(bl_end);
            }
            ++num;
        }

        /**
         * removes specified element at pos.
         * removes the element at pos.
         * returns an iterator pointing to the following element, if pos pointing to the last element, end() will be returned.
         * throw if the container is empty, the iterator is invalid or it points to a wrong place.
         */
        iterator erase(iterator pos) {
            Block bl;
            getblock(0, bl);
            Address block_pos = 0;
            while (!bl.num || bl.nodes[bl.num - 1] < erase_node) {
                if (bl.next == -1) break;
                assert(bl.num);
                block_pos = bl.next;
                getblock(bl.next, bl);
            }
            int erase_pos = lower_bound(bl.nodes, erase_node);
            Node temp = erase_node;
            delete erase_node;
            Node record;
            for (int i = erase_pos + 1; i < bl.num; ++i) {
                bl.nodes[i - 1] = bl.nodes[i];
            }
            --bl.num;
            if (bl.next != -1) {//merge
                Block bl_next;
                getblock(bl.next, bl_next);
                if (bl.num + bl_next.num <= Nmax) {
                    last_block_address = bl.next;
                    for (int i = 0; i < bl_next.num; ++i) {
                        bl.nodes[i + bl.num] = bl_next.nodes[i];
                    }
                    bl.num += bl_next.num;
                    bl.next = bl_next.next;
                }
            }
            putblock(block_pos, bl);
            --num;        }

        /**
         * adds an element to the end
         */
        void push_back(const T &value) {
            insert(end(), value);
        }

        /**
         * removes the last element
         *     throw when the container is empty.
         */
        void pop_back() {
            if (empty())
                throw container_is_empty();
            erase(--end());
        }

        /**
         * inserts an element to the beginning.
         */
        void push_front(const T &value) {
            insert(begin(), value);
        }

        /**
         * removes the first element.
         *     throw when the container is empty.
         */
        void pop_front() {
            erase(begin());
        }
    };

}

#endif
