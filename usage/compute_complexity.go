package main

import (
    "fmt"
    "os"
    comp "github.com/vtphan/repeat-complexity"
)

func main(){
   if len(os.Args) < 2 {
     panic("must provide sequence file.")
   }
   idx := new(comp.Index)
   idx.Build(os.Args[1:])

   k:=6
   fmt.Printf("%s\tR%d\n", os.Args[1], k)
   rks := idx.Rk(k)
   for i:=0; i<len(rks); i++ {
      fmt.Println(i, rks[i])
   }
   idx.Print()
}