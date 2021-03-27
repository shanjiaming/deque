#ifndef SJTU_DEQUE_HPP
#define SJTU_DEQUE_HPP

#include "exceptions.hpp"

#include <cstddef>

#include <cassert>//FIXME delete me
#include <vector>//FIXME delete me
using std::vector;//FIXME delete me

namespace sjtu {

    template<typename T>
    class vector {
    private:
        T** a;
//T**a;//FIXME MemLeak
        size_t space, num;
        static const size_t InitSize = 5;

        void doubleSpace() {
            T** tmp;
//T**tmp;//FIXME MemLeak
            space *= 2;
            tmp = new T*[space];
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
        class const_iterator;
        class iterator {
            friend vector<T>;
        private:
            const vector<T> *vec;
            int pos;
        public:
            iterator(const vector<T> *v, const int p) : vec(v), pos(p) {}

            /*
             * */
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
            iterator& operator+=(const int &n) {
                pos = pos + n;
                return *this;
            }
            iterator& operator-=(const int &n) {
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
            iterator& operator++() {
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
            iterator& operator--() {
                return *this -= 1;
            }
            /**
             * TODO *it
             */
            T& operator*() const{
                return *(vec->a[pos]);
            }
            /**
             * a operator to check whether two iterators are same (pointing to the same memory address).
             */
            bool operator==(const iterator &rhs) const {
                return rhs.vec == vec && rhs.pos == pos;
            }
            bool operator==(const const_iterator &rhs) const {
                return rhs.vec == vec && rhs.pos == pos;
            }
            /**
             * some other operator for iterator.
             */
            bool operator!=(const iterator &rhs) const {
                return !(rhs == *this);
            }
            bool operator!=(const const_iterator &rhs) const {
                return !(rhs == *this);
            }
        };
        /**
         * TODO
         * has same function as iterator, just for a const object.
         */
        class const_iterator {
            friend vector<T>;
        private:
            const vector<T> *vec;
            int pos;
        public:
            const_iterator(const vector<T> *v, const int p) : vec(v), pos(p) {}

            /**
             * return a new const_iterator which pointer n-next elements
             * as well as operator-
             */
            const_iterator operator+(const int &n) const {
                return const_iterator(vec, pos + n);
            }
            const_iterator operator-(const int &n) const {
                return const_iterator(vec, pos - n);
            }
            // return the distance between two iterators,
            // if these two iterators point to different vectors, throw invaild_iterator.
            int operator-(const const_iterator &rhs) const {
                if (vec != rhs.vec)
                    throw invalid_iterator();
                return pos - rhs.pos;
            }
            const_iterator& operator+=(const int &n) {
                pos = pos + n;
                return *this;
            }
            const_iterator& operator-=(const int &n) {
                pos = pos - n;
                return *this;
            }
            /**
             * TODO iter++
             */
            const_iterator operator++(int) {
                auto tmp = *this;
                *this += 1;
                return tmp;
            }
            /**
             * TODO ++iter
             */
            const_iterator& operator++() {
                return *this += 1;
            }
            /**
             * TODO iter--
             */
            const_iterator operator--(int) {
                auto tmp = *this;
                *this -= 1;
                return tmp;
            }
            /**
             * TODO --iter
             */
            const_iterator& operator--() {
                return *this -= 1;
            }
            /**
             * TODO *it
             */
            const T& operator*() const{
                return *(vec->a[pos]);
            }
            /**
             * a operator to check whether two iterators are same (pointing to the same memory address).
             */
            bool operator==(const iterator &rhs) const {
                return rhs.vec == vec && rhs.pos == pos;
            }
            bool operator==(const const_iterator &rhs) const {
                return rhs.vec == vec && rhs.pos == pos;
            }
            /**
             * some other operator for const_iterator.
             */
            bool operator!=(const iterator &rhs) const {
                return !(rhs == *this);
            }
            bool operator!=(const const_iterator &rhs) const {
                return !(rhs == *this);
            }
        };
        /**
         * TODO Constructs
         * Atleast two: default constructor, copy constructor
         */
        vector() : space(InitSize), num(0) {
            a = new T*[InitSize];
        }
        vector(const vector &other) : space(other.space), num(other.num) {
            a = new T*[space];
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
            space = other.space;
            num = other.num;
            a = new T*[space];
            for (int i = 0; i < num; ++i) {
                a[i] = new T(*other.a[i]);
            }

            return *this;
        }
        /**
         * assigns specified element with bounds checking
         * throw index_out_of_bound if pos is not in [0, space)
         */
        T & at(const size_t &pos) {
            if (pos < 0 || pos >= num)
                throw index_out_of_bound();
            return *a[pos];
        }
        const T & at(const size_t &pos) const {
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
        T & operator[](const size_t &pos) {
            if (pos < 0 || pos >= num)
                throw index_out_of_bound();
            return *a[pos];
        }
        const T & operator[](const size_t &pos) const {
            if (pos < 0 || pos >= num)
                throw index_out_of_bound();
            return *a[pos];
        }
        /**
         * access the first element.
         * throw container_is_empty if space == 0
         */
        const T & front() const {
            if (!num)
                throw container_is_empty();
            return *a[0];
        }
        /**
         * access the last element.
         * throw container_is_empty if space == 0
         */
        const T & back() const {
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
            return insert(pos.pos, value);
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
            return erase(pos.pos);
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

        using Node = T;

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
            b = blocks[x];
        }

        void putblock(Address x, const Block &b) {
            blocks[x] = b;
        }

        void putblockend(const Block &b) {
            blocks.push_back(b);
        }


    public:
        class iterator;

        typedef iterator const_iterator;

//	class const_iterator;
        class iterator {
            friend deque<T>;
        public:

            deque *deque_ptr;
            int block_th;
            int node_th;
            int endnote;
            int dist;

            iterator() {};

            iterator(deque *_deque_ptr, int _block_th, int _node_th, int _endnote, int _dist) : deque_ptr(_deque_ptr),
                                                                                         block_th(_block_th),
                                                                                         node_th(_node_th),
                                                                                         endnote(_endnote) ,
                                                                                         dist(_dist){
            };

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
                if (deque_ptr != rhs.deque_ptr) throw invalid_iterator();
                return dist - rhs.dist;
            }

            iterator &operator+=(const int &n) {
                //TODO
                if (n < 0)
                    return (*this -= (-n));
                dist += n;
                int nchange = n + node_th + endnote;
                endnote = 0;
                while (block_th < deque_ptr->blocks.size()) {
                    if (nchange < deque_ptr->blocks[block_th].num) {
                        node_th = nchange;
                        return *this;
                    }
                    nchange -= deque_ptr->blocks[block_th].num;
                    ++block_th;
                }
                endnote = nchange;
                node_th = 0;
                return *this;
            }

            iterator &operator-=(const int &n) {
                if (n < 0)
                    return (*this += (-n));
                dist -= n;
                int nchange = -n + node_th + endnote;
                endnote = 0;
                while (block_th >= 0) {
                    if (nchange >= 0) {
                        node_th = nchange;
                        return *this;
                    }
                    --block_th;
                    nchange += deque_ptr->blocks[block_th].num;
                }
                endnote = nchange;
                node_th = 0;
                return *this;
            }


            iterator operator++(int) {
                iterator ret = *this;
                *this += 1;
                return ret;
            }

            iterator &operator++() {
                return *this += 1;
            }

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
                if (endnote != 0) throw invalid_iterator();
                return (deque_ptr->blocks[block_th].nodes[node_th]);
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

//        deque(const deque &other) {
//            blocks = other.blocks;
//            last_block_address = other.last_block_address;
//            num = other.num;
//        }

        deque(const deque &other) = default;

        /**
         * TODO Deconstructor
         */
        ~deque() {}

        /**
         * TODO assignment operator
         */
//        deque &operator=(const deque &other) {
//            blocks = other.blocks;
//            last_block_address = other.last_block_address;
//            num = other.num;
//            return *this;
//        }

        deque &operator=(const deque &other) = default;

        /**
         * access specified element with bounds checking
         * throw index_out_of_bound if out of bound.
         */
        T &at(const size_t &pos) {
            try {
                auto& ret = *(begin() + pos);
                return ret;
            } catch (invalid_iterator) {
                throw index_out_of_bound();
            }
        }

        const T &at(const size_t &pos) const {
            try {
                auto ret = *(cbegin() + pos);
                return ret;
            } catch (invalid_iterator) {
                throw index_out_of_bound();
            }
            //return *(cbegin() + pos);//FIXME ����汾��ʱ������ôд����Ϊconst_iterator�Ǽ�const��
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
            return iterator(this, 0, 0, 0, 0);
        }

        const_iterator cbegin() const {}

        /**
         * returns an iterator to the end.
         */
        iterator end() {
            return iterator(this, last_block_address, blocks[last_block_address].num - 1, 1, num);
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
            if (pos.endnote == 1) {
                push_back(value);
                return --end();
            }
            if (pos.endnote != 0) {
                throw invalid_iterator();
            }
            Address block_pos = pos.block_th;
            Block &bl = blocks[block_pos];
            int insert_pos = pos.node_th;

            if (bl.num < Nmax) {//insert directly
                if (bl.num == insert_pos) {
                    bl.nodes.push_back(value);
//                bl.nodes[insert_pos] = value;
                } else {
                    Node temp = value;
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
                    bl.nodes[insert_pos] = value;
                    bl_end.num = bl.num - insert_pos;
                    bl.num = insert_pos + 1;
                } else {
                    bl_end.nodes.push_back(value);
                    bl_end.num = bl.num - insert_pos + (insert_pos == bl.num);
                }
                bl_end.next = bl.next;
                bl.next = blocks.size() - 1;
                if (bl_end.next == -1) {
                    last_block_address = bl.next;
                }
                //this is change mark
                putblock(block_pos, bl);
                putblockend(bl_end);
            }
            ++num;
            return pos;//FIXME
        }

        /**
         * removes specified element at pos.
         * removes the element at pos.
         * returns an iterator pointing to the following element, if pos pointing to the last element, end() will be returned.
         * throw if the container is empty, the iterator is invalid or it points to a wrong place.
         */
        iterator erase(iterator pos) {
            if (pos + 1 == end()) {
                pop_back();
                return end();
            }
            if (pos.endnote != 0)
                throw invalid_iterator();
            if (empty())
                throw container_is_empty();
            Address block_pos = pos.block_th;
            Block &bl = blocks[block_pos];
            int erase_pos = pos.node_th;
            Node record;
            for (int i = erase_pos + 1; i < bl.num; ++i) {
                bl.nodes[i - 1] = bl.nodes[i];
            }
            --bl.num;
            if (bl.next != -1) {//merge
                Block bl_next;
                getblock(bl.next, bl_next);
                if (bl.num + bl_next.num <= Nmax) {
                    if (bl_next.next == -1)
                        last_block_address = bl.next;//FIXME endϵ�е��ж���charu������ĳ�pushbsack��Ҳ�������⡣���й���redundancy�ļ��㣨endnote��
                    for (int i = 0; i < bl_next.num; ++i) {
                        bl.nodes[i + bl.num] = bl_next.nodes[i];
                    }
                    bl.num += bl_next.num;
                    bl.next = bl_next.next;
                }
            }
            putblock(block_pos, bl);
            --num;
            return pos;//FIXME
        }

        /**
         * adds an element to the end
         */
        void push_back(const T &value) {
            if (last_block_address == -1){
                blocks.push_back(Block());
                last_block_address = blocks.size()  - 1;
            }
            Block &bl = blocks[last_block_address];
            assert(bl.next == -1);
            ++num;
            if (bl.num < Nmax) {
                bl.nodes.push_back(value);
                bl.num++;
                return;
            }
            Block bl_end;
            bl_end.nodes.push_back(value);
            bl_end.num = 1;
            bl_end.next = -1;
            blocks.push_back(bl_end);
            last_block_address = bl.next = blocks.size() - 1;
        }

        /**
         * removes the last element
         *     throw when the container is empty.
         */
        void pop_back() {
            if (empty())
                throw container_is_empty();
            --num;
            Block &bl = blocks[last_block_address];
            assert(bl.next == -1);
            bl.nodes.pop_back();
            bl.num--;
            if (bl.num > 1) {
                return;
            }
            //kill the last block, O(n)
            Address bl_pos = 0, bl_next_pos = blocks[bl_pos].next;
            while(blocks[bl_next_pos].next != -1){
                bl_pos = bl_next_pos;
                bl_next_pos = blocks[bl_next_pos].next;
            }
            last_block_address = bl_pos;
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
