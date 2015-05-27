/*
Author: Vinhthuy Phan, Shanshan Gao
Copyright 2014
Measures of complexity: I, Ik, D, Dk, Rk
*/
package genomecomplexity

import (
    // "fmt"
    "os"
    "bufio"
    "bytes"
    "io/ioutil"
    "fmt"
)

type Index struct{
    data []byte
    sa []int
    lcp []int
    Length int
    acc_length []int      // accumulated lengths of consecutive contigs
    Repeat_count []int    // number of global repeats in each contig
}

func (idx *Index) Print() {
  // fmt.Println(string(idx.data), idx.Length, idx.acc_length)
  // fmt.Println(idx.sa)
  // fmt.Println(idx.lcp)
  fmt.Println(idx.Repeat_count)
  fmt.Println(idx.acc_length)
}

func (idx *Index) Build(filenames []string) {
  idx.acc_length = make([]int, len(filenames))
  idx.Repeat_count = make([]int, len(filenames))
  for j:=0; j<len(filenames); j++ {
    if j==0 {
      idx.data = ReadSequence(filenames[0])
    } else {
      idx.data = append(idx.data, 'N')
      idx.data = append(idx.data, ReadSequence(filenames[j])...)
    }
    idx.acc_length[j] = len(idx.data)
  }
  ws := &WorkSpace{}
  idx.sa = make([]int, len(idx.data))
  ws.ComputeSuffixArray(idx.data, idx.sa)
  // idx.sa = qsufsort(idx.data)
  idx.lcp = make([]int, len(idx.data)-1)  // lcp[i] stores length of lcp of sa[i] and sa[i+1]
  for i := 1; i < len(idx.data); i++ {
    idx.lcp[i-1] = idx.lcp_len(i)
  }
  idx.Length = len(idx.sa)
}

// length of longest common prefix of data[SA[m]:] and data[SA[m-1]:]
func (idx *Index) lcp_len(m int) int{
    L, i, j := len(idx.data), idx.sa[m], idx.sa[m-1]
    for i<L && j<L && idx.data[i]==idx.data[j] {
        i++
        j++
    }
    return j - idx.sa[m-1]
}

// return the interval containing a k-mer starting at position i ([i, i+k-1])
func (idx Index) Locate_contig(i int, k int) int {
  begin := 0
  for j:=0; j<len(idx.acc_length); j++ {
    if i<begin {
      return -1
    }
    if i<=idx.acc_length[j]-k {
      return j
    }
    begin = idx.acc_length[j]+1
  }
  return -1
}

func (idx Index) Block(m int, k int) int {
  for i := m; i < len(idx.data)-1; i++ {
    fmt.Println("\t",i, idx.sa[i], string(idx.data[idx.sa[i] : idx.sa[i]+k]), "contig", idx.Locate_contig(idx.sa[i], k))
    if idx.lcp[i] < k {
        return (i - 1)
    }
  }
  return len(idx.data) - 2
}

// Rk = k-repeat density
func (idx Index) Rk(k int) []float64{
  var c uint64 = 0
  i, j := 0, 0
  var contig int
  contig_flag := make([]bool, len(idx.Repeat_count))

  for i < len(idx.data)-1 {
    // fmt.Println(i, idx.lcp[i], idx.Block(i,k), c)
    if idx.lcp[i] >= k {
      j = idx.Block(i,k)

      for x:=0; x<len(contig_flag); x++ {
        contig_flag[x] = false
      }
      for x:=i; x<=j+1; x++ {
        contig = idx.Locate_contig(idx.sa[x], k)
        if contig >= 0 {
          contig_flag[contig] = true
        }
      }
      for x:=0; x<len(contig_flag); x++ {
        if contig_flag[x] {
          idx.Repeat_count[x] += (j-i+2)
        }
      }

      c += uint64(j - i + 2)
      i = j + 1
    } else {
      i++
    }
  }
  contig_rk := make([]float64, len(idx.Repeat_count))
  prev := 0
  for i:=0; i<len(contig_rk); i++ {
    contig_rk[i] = float64(idx.Repeat_count[i])/float64(idx.acc_length[i] - prev - k + 1)
    prev = idx.acc_length[i]+1
  }
  return contig_rk
}

func ReadSequence(file string) []byte{
   f, err := os.Open(file)
   if err != nil {
      panic(err)
   }
   defer f.Close()
   byte_array := make([]byte, 0)
   Ns := []byte("N")
   None := []byte("")
   if file[len(file)-6:] == ".fasta" {
      scanner := bufio.NewScanner(f)
      for scanner.Scan() {
         line := scanner.Bytes()
         if len(line)>0 && line[0] != '>' {
            byte_array = append(byte_array, bytes.Replace(bytes.Trim(line,"\n\r "), Ns, None, -1)...)
         }
      }
   } else {
      byte_array, err = ioutil.ReadFile(file)
      if err != nil {
         panic(err)
      }
   }
   return byte_array
}