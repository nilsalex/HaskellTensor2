{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE LambdaCase #-}

module LinearAlgebraFF (linearAlgebraFFTest) where

import Test.QuickCheck
import Test.Tasty.QuickCheck

import Test.Tasty

import System.Exit

import Numeric.LinearAlgebra (rank)
import qualified Numeric.LinearAlgebra.Data as Matrix

import Math.Tensor.Internal.LinearAlgebra (independentColumnsFF, independentColumnsMatFF, gaussianFF, isref, verify)

data SmallInt = S0 | S1 deriving (Show, Ord, Eq, Enum, Bounded)

rankFF :: Matrix.Matrix Matrix.Z -> Int
rankFF mat = rank mat'
  where
    mat' = Matrix.fromZ mat :: Matrix.Matrix Double

toSmall :: Int -> SmallInt
toSmall 0 = S0
toSmall 1 = S1
toSmall i = error $ "cannot convert " ++ show i ++ " to SmallInt"

fromSmall :: Num a => SmallInt -> a
fromSmall S0 = 0
fromSmall S1 = 1

instance Arbitrary SmallInt where
    arbitrary = arbitraryBoundedEnum

data MatrixData a = MatrixData (Positive Int) (Positive Int) [a] deriving Show

instance Arbitrary a => Arbitrary (MatrixData a) where
    arbitrary = do
      m@(Positive m') <- arbitrary
      n@(Positive n') <- arbitrary
      xs <- vector (m'*n')
      return $ MatrixData m n xs

prop_smallValues :: MatrixData SmallInt -> Bool
prop_smallValues (MatrixData (Positive rows) (Positive cols) xs) =
    rankFF mat' == rankFF mat &&
    verify mat ref
  where
    mat  = (rows Matrix.>< cols) $ map fromSmall xs
    mat' = independentColumnsMatFF mat
    ref  = gaussianFF mat

prop_ints :: MatrixData Int -> Bool
prop_ints (MatrixData (Positive rows) (Positive cols) xs) =
    rankFF mat' == rankFF mat &&
    verify mat ref
  where
    mat  = (rows Matrix.>< cols) $ map fromIntegral xs
    mat' = independentColumnsMatFF mat
    ref  = gaussianFF mat

prop_consec :: Positive Int -> Int -> Bool
prop_consec (Positive dim') start =
    icols == [0,1] &&
    verify mat ref
  where
    dim   = dim' + 100
    mat   = (dim Matrix.>< dim) $ map fromIntegral [start..]
    icols = independentColumnsFF mat
    ref   = gaussianFF mat

testCase1 = localOption (QuickCheckMaxSize 10) $ localOption (QuickCheckTests 100000) $ testProperty "prop_smallValues" prop_smallValues
testCase2 = localOption (QuickCheckMaxSize 10) $ localOption (QuickCheckTests 100000) $ testProperty "prop_ints" prop_ints
testCase3 = localOption (QuickCheckTests 1000) $ testProperty "prop_consec" prop_consec

linearAlgebraFFTest = testGroup "LinearAlgebraFFTest" [testCase1, testCase2, testCase3]
