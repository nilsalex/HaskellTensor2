module Main (main) where

import Ansatz
import IndList
import LinearAlgebra
import LinearAlgebraFF
import Serialization

import Test.Tasty

test = testGroup "sparse-tensor tests"
        [
          indListTest,
          linearAlgebraTest,
          linearAlgebraFFTest,
          ansatzTest,
          serializationTest
        ]

main :: IO ()
main = defaultMain test
