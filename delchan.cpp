/* 
 * File:   main.cpp
 * Author: Jarek Duda
 * Created on July 27, 2014, 12:31 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include<vector>
#include<algorithm>

using namespace std;
typedef unsigned char byte;

const short R=7, Rc = 8-R, Rcp = 1 << (Rc + (R==7));        // R is the number of redundancy bits per byte (rate = Rc/8)
const uint64_t rmask = (1<<R)-1, xmask = (1<<Rc)-1; 
const uint64_t INITSTATE = 0;                       // initial state for the coding
const int STEP_LIMIT = 5E7;                         // the main tree will use STEP_LIMIT * 23 bytes

int count1(byte c) {int s=0; while(c){s+= c&1; c>>=1;}; return s;}   // count "1"
void printc(uint64_t c, int n=8) {for(int i=0;i<n;i++) cout<< (1&(c>>i)); cout<<" ";}

struct coding{                      // THE CODING CONTAINER   
    uint64_t tr[Rcp];               // transition function - defines encoding
    uint64_t state;                 // current state    
    coding(short R){                // generate pseudorandom tr[] defining the coding                
        for(int i=0 ; i < Rcp ; i++)                     // block code for the first block
	{   if (R == 1) tr[i] = 1&((i>>6)^(i>>5)^(i>>4)^(i>>3)^(i>>2)^(i>>1)^i);
            if (R == 2) tr[i] = 2*(1&((i>>3)^(i>>2)^(i>>1)^i)) + (1&((i>>5)^(i>>4)^(i>>3)^(i>>2)));
            if (R == 3) tr[i] = 4*(1&((i>>3)^(i>>2)^(i>>1)^i)) + 2*(1&((i>>4)^(i>>2)^(i>>1)^i)) + (1&((i>>4)^(i>>3)^(i>>1)^i));
            if (R == 4) tr[i] = 8*(1&((i>>2)^(i>>1)^i)) + 4*(1&((i>>3)^(i>>1)^i)) + 2*(1&((i>>3)^(i>>2)^i)) + (1&((i>>3)^(i>>2)^(i>>1)));
            if (R == 5) tr[i] = (1&(i^(i>>1)))+(2&(i^(i>>1)))+4*(1&(i^(i>>2)))+24*(1&(i^(i>>1)^(i>>2)));
            if (R == 6) tr[i] = 3*(1&i)+6*(2&i)+48*(1&(i^(i>>1)));
            if (R == 7) tr[i] = 127*i;
	}
        uint64_t m = 1<<R;
	int temp[Rcp];
        for(int i = 0 ; i < (64-R)/ Rc ; i++)              // fill the remaining bits
	{   for (int j=0 ; j < Rcp ; j++) temp[j] = j;        // find a random permutation
	    for (int32_t j=0 ; j < Rcp ; j++)
	    {  	int k = floor((float)rand()/RAND_MAX *(Rcp - j));
	    	tr[j] += m*temp[k];
                temp[k] = temp[Rcp - j - 1];  
	     }  m<<=Rc;        
        }
    }
    inline byte encstep(byte symbol){              // encoding step
        uint64_t t=tr[symbol >> R];
        byte y = symbol | ((state ^ t) & rmask);
        state = ((state ^ t) >> R) | (state << 64-R);
        return y;
    }
    inline bool decstep(byte symbol){              // decoding step
        uint64_t t=tr[symbol >> R];
        bool correct = (((symbol ^ state ^ t) & rmask) == 0);            // if redundancy bits are correct
        state = ((state ^ t) >> R) | (state << (64-R));
        return correct;
    }
    vector<byte> encode(vector<byte> & msg, uint64_t & last_state){         // encode message
        vector<byte> enc;
        uint32_t buf=0, nb=0; state=INITSTATE;
            for(byte c:msg){                              // encode
                buf += c<<nb; nb+=8;                // read byte to buffer
                while(nb >= Rc){                    
                    enc.push_back(encstep((buf & xmask) << R));             
                    buf >>= Rc; nb -= Rc;
                }
            }
            if(nb) enc.push_back(encstep(buf<<R));  // empty buffer
            last_state = state;
            return enc;
    }
    vector<byte> decode(vector<byte> & enc){        // pure decode (no correction))
        vector<byte> dec;
        uint32_t buf=0, nb=0; state=INITSTATE;
        for(byte c:enc){
            buf += (c >> R) << nb; nb += Rc;            
            if(nb>=8){dec.push_back(buf & 255); buf >>= 8; nb-=8; }            
        }
        return dec;
    }
};

class test{                             // GENERATING MESSAGE, ENCODING, DAMAGING AND THEN CHECKING
    int L;  // length
    vector<byte> msg, enc, dec;   // message, encoded message, damage vector, damaged encoded message
public:
    vector<byte> dam;
    int len;                              // number of bits in damaged message
    uint64_t last_state;
    test(coding & cod, int ln, float pr){              // constructor - generate message, encode, damage
        L=ln; 
        for(int i=0;i<L;i++) msg.push_back(255&rand());       //generate random message
        enc = cod.encode(msg, last_state);                    // encode
        uint32_t buf=0, nb=0;                                // damage enc
        for(int i=0; i<enc.size(); i++){
            byte c=enc[i];
            for(int j=0;j<8;j++){
                if((float)rand()/RAND_MAX > pr) buf += (c&1) << (nb++);  // copy if not deleted
                c >>= 1;
            }
            if(nb>8) {dam.push_back(buf&255); buf >>= 8; nb -= 8;};
        }
        len = 8*dam.size() + nb;           // number of bits in damaged message
        if(nb) dam.push_back(buf);         // eventually empty buffer
        //for(int i=0; i<enc.size(); i++)  {printc(enc[i]); cout<<"("<<(int)enc[i]<<"), dam: "; printc(dam[i]); cout<<" ("<<(int)dam[i]<<")"<<endl;}
    }
    float check(vector<byte> & rep){                    // what fraction of bytes was properly corrected
        float agree = 0;
        for(int i=0; i<min(rep.size(), enc.size()); i++) if(rep[i]==enc[i]) agree++;
        return agree/enc.size();
    }
};

// -------------------- THE CORRECTION PART -------------------------------
 struct corlst {            // for ~1MB list of all corrections (at least 1 bit deleted)
     byte symbol;           // corrected symbol;
     byte nb;               // how many bits of the stream it used (7,6,5,4,3,2,1,0)
     float weight;           // weight of using this correction
     bool operator<(const struct corlst & b) const { return this->weight > b.weight;}
 };
struct Node {               // node of the correction tree
  uint16_t position;        // currently decoded bit block - change to 32 if larger than 8kB frame         
  uint64_t state;           // current state
  uint32_t parent;          // the number of parent node
  float weight;             // weight of this node
  uint32_t last_child;      // last considered child for this node
  byte symbol;              // symbol that lead the parent here  
};
struct HeapNode {                   // we need heap of nodes worth considering
  uint32_t node;                    // which node it corresponds to
  float weight;                     // weight after applying nodes[node].last_child 
  bool operator<(const struct HeapNode & b) const { return this->weight < b.weight;}
};
typedef vector<HeapNode> Heap;
inline void heap_push(Heap & heap, HeapNode n){
    heap.push_back(n); push_heap(heap.begin(),heap.end()); }
inline HeapNode heap_pop(Heap & heap){
    HeapNode n = heap.front();pop_heap(heap.begin(),heap.end());heap.pop_back(); return n;
}
struct BSTnode{
    uint64_t val;
    uint32_t left, right;
}; 
class BST{                              // to test if haven't already reached this position with the same state
    vector<BSTnode> tree;
public:
    bool spotted(uint64_t v){           // tries to add to BST and returns true if already spotted there
        if(tree.empty()) {tree.push_back({v,0,0}); return false;} else{
            uint32_t cn=0, pn;
            bool found = true;
            while(v != tree[cn].val){
                pn = cn;
                if(v < tree[cn].val){
                    cn=tree[cn].left; 
                    if(cn==0) {found = false; cn = tree[pn].left=tree.size(); tree.push_back({v,0,0});}
                }
                else {
                    cn=tree[cn].right; 
                    if(cn==0) {found = false; cn = tree[pn].right=tree.size(); tree.push_back({v,0,0});}
                }
            }
        return found;
        }        
    }
};

class correct{                                    // tables and reservations to speedup tests 
    vector<corlst> error_list = {{0,0,0}};        // lists of all corrections
    uint32_t error_start[128][1<<R];              // for all blocks and last R bits of the state
    vector<HeapNode> heap;                        // heap of nodes worth considering
    vector<Node> nodes;                           // huge table of all nodes
    float pr;
public:
    int steps;
    correct(coding &cod, float p){
        pr=p; heap.reserve(1 + (STEP_LIMIT >> 5));nodes.reserve(1 + (STEP_LIMIT >> 5));
        float prb[9]; for(int nb=0; nb<9; nb++) prb[nb] = log(pow(pr,nb) * pow(1-pr,8-nb));
        vector<byte> tnb[9];                      // tnb[i] is all bytes with i of "1"     
        for(int i=0; i<256; i++) tnb[count1(i)].push_back(i);
        vector<byte> cor;                           // current correction list
        vector<pair<byte,uint16_t>> crl[1<<R];      // list of corrections    
   
        // generating error_list and error_start[b][red] for current block b and first state bits red
        for(byte b=0; b<128; b++){
            for(int i=0; i < (1<<R); i++) crl[i].clear();
            for(int nb=1; nb<9; nb++){
                cor.clear();                    // all possible corrections for block b and nb deletions
                for(int msk:tnb[nb]){
                    for(int j=0; j < (1<<nb); j++){
                        int pb=7-nb, pj=nb-1, fin=0;
                        for(int k=0; k < 8; k++){ 
                            if(msk & (1 << k)) {fin = 2*fin + (1 & (j >> pj)); pj--;} else 
                                {fin = 2*fin +  (1 & (b >> pb)); pb --;}             
                        } cor.push_back(fin); 
                    }
                }
                sort(cor.begin(),cor.end());        // to count duplicating corrections
                cor.push_back(cor[0]);              //guardian
                int cc=cor[0], cn=1;                // cn is the number of appearances of given correction
                for(int i=1; i<cor.size(); i++){
                    if(cor[i]==cc) cn++; else {                        
                        crl[(cod.tr[cc>>R] ^ cc) & rmask].push_back({cc, cn + (nb<<8)});                    
                        cc=cor[i]; cn=1;
                    }
                }
            }               
            for(int i=0; i < (1<<R); i++){
                error_start[b][i]=error_list.size();
                for(auto a:crl[i]){
                    byte nb = a.second >> 8, cn = (a.second & 255);      // decode
                    error_list.push_back({a.first, 8-nb, prb[nb] + log(2)*R });  // + log(cn)? no!
                }
                // sort(error_list.begin() + error_start[b][i], error_list.end());   // not needed without cn
                error_list.push_back({b,0,-100000});                              // guardian
            }
        }
    }
    vector<byte> repair(coding & cod, vector<byte> dam, int len, uint64_t last_state){  // the repair procedure
        byte blocks[len+1];           // blocks contains all succeeding 8 bit blocks (overlapping))
        uint32_t buf=0, nb=0, i=0; 
        dam.push_back(0);               // for eventual deletions at the end 
        for(byte c:dam){                              
            buf += c<<nb; nb+=8;
            while(nb>8) {blocks[i++] = (buf & 255); buf >>=1; nb --;}
        }                                                  
        heap.clear(); nodes.clear();
        vector<BST> bst(len+1);                               // table of roots of BST tree for all positions
        const float wzero = 8*log(1-pr) + log(2)*R;           // weight for node without error    
        Node nn;                                              // new node
        nn.last_child=1;    nodes.push_back(nn);              // root and guardian
        nn.last_child=0; 
        volatile bool correct = false;
        cod.state = INITSTATE;                                // preparing the root of the correction tree
        uint32_t parent = 0, position = 0, cpos, ppos;                                                     
        byte symbol = 0;                                                               
        float weight = 0;                                      
      // --------------------  THE CORRECTION LOOP -------------------------
      while (!correct){                                                         
          correct = true;                               // there are only nodes passing verification here   
          while (correct && position <= len && nodes.size() < STEP_LIMIT ) {       // create node and try to expand until error
             nn.position = position;                                 
             nn.symbol = symbol;                                    // parent's symbol leading to this node
             nn.state = cod.state;                                  // state in this node
             nn.parent = parent;                                    // parent is nodes[parent]
             nn.weight = weight;                                    // weight here
             symbol = blocks[position];      // will try error-free branch (the most probable)          
             nn.last_child = error_start[symbol & 127][nn.state & rmask];
             parent = nodes.size();             
             nodes.push_back(nn);                                    // store the new node in the node table
             weight += wzero;                // Increase the weight by a error-free weight
             correct = cod.decstep(symbol);  // Check if the uncorrected symbol passes verification
             heap_push(heap, {parent, nn.weight + error_list[nn.last_child].weight});      // insert first correction to heap
             position += 8;                  // assuming no delation   
         }
         if(position > len) correct = false;           // getting out of message with wrong state
         if(nn.state == last_state) correct = true;             // successful correction
         if(nodes.size() >= STEP_LIMIT) correct = true;         // emergency exit if out of memory
         if(!correct)                       // the loop exit was due to an error
             do{HeapNode tn = heap_pop(heap);               // get new node
                parent = tn.node;
                Node * temp = & nodes[parent];
                symbol = error_list[temp->last_child].symbol;
                position = temp->position + error_list[temp->last_child].nb;
                heap_push(heap, {parent, temp->weight + error_list[++ (temp->last_child)].weight});   // push next correction        
                weight = tn.weight;
                cod.state = temp->state;
                cod.decstep(symbol);            
            } while(position <= len && bst[position].spotted(cod.state));              // check if this state hasn't been already reached here
    }
     vector<byte> rep ;                      // apply the correction 
     while(parent > 1){nn = nodes[parent]; parent = nn.parent; rep.push_back(nn.symbol);}
     reverse(rep.begin(), rep.end());
     steps = nodes.size();
     heap.clear(); nodes.clear();
     return rep;   
    }        
};

int main() {   
    ofstream out("dataR7");                     // R is chosen at the top (rate = (8-R)/8)
    srand(123); //srand(time(0));
    coding cod(R);
    for(float pr=0.1; pr<0.5; pr+=0.01){       // range of deletion probabilities to test
       correct cor(cod,pr);                     // generate correction tables
       for(int j=0; j<1000; j++){               // 1000 tests per probability
           srand(j);
           test t(cod, 125 * Rc, pr);           // 1000 byte encoded message
           vector<byte> repaired = cor.repair(cod, t.dam, t.len, t.last_state); 
           out<<pr<<" "<<j<<" "<<cor.steps<<" "<<t.check(repaired)<<endl;    
           //cout<<pr<<" "<<j<<" "<<cor.steps<<" "<<t.check(repaired)<<endl; 
       }   
    }        
    out.close();
    return 0;
}
