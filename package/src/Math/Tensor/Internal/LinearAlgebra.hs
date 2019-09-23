-----------------------------------------------------------------------------
-- |
-- Module      :  Math.Tensor.Internal.LinearAlgebra
-- Copyright   :  (c) 2019 Tobias Reinhart and Nils Alex
-- License     :  MIT
-- Maintainer  :  tobi.reinhart@fau.de, nils.alex@fau.de
--
-- Gaussian elimination algorithm based on hmatrix.
-----------------------------------------------------------------------------

{-# LANGUAGE ScopedTypeVariables #-}

module Math.Tensor.Internal.LinearAlgebra (
-- * Gaussian Elimination
gaussianST, gaussianFFST,
gaussian, gaussianFF,
-- * Linearly Independent Columns
independentColumns, independentColumnsFF,
independentColumnsMat, independentColumnsMatFF,
-- * Pivots
pivotsU, pivotsUFF,
findPivotMax, findPivotMaxFF)

where

import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Numeric.LinearAlgebra.Devel

import Data.List (maximumBy)

import Control.Monad
import Control.Monad.ST

-- | Returns the pivot columns of an upper triangular matrix.
--
-- @
-- &#x3BB; let mat = (3 >< 4) [1, 0, 2, -3, 0, 0, 1, 0, 0, 0, 0, 0]
-- &#x3BB; mat
-- (3><4)
--  [ 1.0, 0.0, 2.0, -3.0
--  , 0.0, 0.0, 1.0,  0.0
--  , 0.0, 0.0, 0.0,  0.0 ]
-- &#x3BB; pivotsU mat
-- [0,2]
-- @
--

pivotsU :: Matrix Double -> [Int]
pivotsU mat = go (0,0)
  where
    go (i,j)
      = case findPivot mat e (i,j) of
          Nothing       -> []
          Just (i', j') -> j' : go (i'+1, j'+1)
    maxAbs = maximum $ map (maximum . map abs) $ toLists mat
    e = eps * maxAbs

pivotsUFF :: Matrix Z -> [Int]
pivotsUFF mat = go (0,0)
  where
    go (i,j)
      = case findPivotFF mat (i,j) of
          Nothing       -> []
          Just (i', j') -> j' : go (i'+1, j'+1)

eps :: Double
eps = 1e-12

-- find next pivot in upper triangular matrix

findPivotFF :: Matrix Z -> (Int, Int) -> Maybe (Int, Int)
findPivotFF mat (i, j)
    | n == j = Nothing
    | m == i = Nothing
    | otherwise = case nonZeros of
                    []          -> if n == j+1
                                   then Nothing
                                   else findPivotFF mat (i, j+1)
                    (pi, pj):_  -> Just (pi, pj+j)
    where
        m = rows mat
        n = cols mat
        col = mat ¿ [j]
        nonZeros = filter (\(i', _) -> i' >= i) $ find (/= 0) col

findPivot :: Matrix Double -> Double -> (Int, Int) -> Maybe (Int, Int)
findPivot mat e (i, j)
    | n == j = Nothing
    | m == i = Nothing
    | otherwise = case nonZeros of
                    []          -> if n == j+1
                                   then Nothing
                                   else findPivot mat e (i, j+1)
                    (pi, pj):_  -> Just (pi, pj+j)
    where
        m = rows mat
        n = cols mat
        col = mat ¿ [j]
        nonZeros = filter (\(i', _) -> i' >= i) $ find (not . (< e) . abs) col

-- | Find pivot element below position (i, j) with greatest absolute value in the ST monad.

findPivotMaxFF :: Int -> Int -> Int -> Int -> STMatrix s Z -> ST s (Maybe (Int, Int))
findPivotMaxFF m n i j mat
    | n == j = return Nothing
    | m == i = return Nothing
    | otherwise =
        do
          col      <- mapM (\i' -> do
                                    x <- readMatrix mat i' j
                                    return (i', abs x))
                      [i..m-1]
          let nonZeros = filter ((/= 0) . abs . snd) col
          case nonZeros of
            []        -> if n == j+1
                         then return Nothing
                         else findPivotMaxFF m n i (j+1) mat
            (pi,_):_  -> return $ Just (pi, j)

findPivotMax :: Int -> Int -> Int -> Int -> STMatrix s Double -> ST s (Maybe (Int, Int))
findPivotMax m n i j mat
    | n == j = return Nothing
    | m == i = return Nothing
    | otherwise =
        do
          col      <- mapM (\i' -> do
                                    x <- readMatrix mat i' j
                                    return (i', abs x))
                      [i..m-1]
          let nonZeros = filter (not . (<eps) . abs . snd) col
          let (pi, _) = maximumBy (\(_, x) (_, y) -> x `compare` y) nonZeros
          case nonZeros of
            [] -> if n == j+1
                  then return Nothing
                  else findPivotMax m n i (j+1) mat
            _  -> return $ Just (pi, j)

-- gaussian elimination of sub matrix below position (i, j)

gaussianFF' :: Int -> Int -> Int -> Int -> STMatrix s Z -> ST s ()
gaussianFF' m n i j mat = do
    iPivot' <- findPivotMaxFF m n i j mat
    case iPivot' of
        Nothing     -> return ()
        Just (r, p) -> do
                          rowOper (SWAP i r (FromCol j)) mat
                          pv <- readMatrix mat i p
                          mapM_ (reduce pv p) [i+1 .. m-1]
                          gaussianFF' m n (i+1) (p+1) mat
  where
    reduce pv p r = do
                      rv <- readMatrix mat r p
                      if rv == 0
                        then return ()
                        else
                         let --frac = -rv / pv
                             op1 = SCAL pv (Row r) (FromCol p)
                             op2 = AXPY (-rv) i r (FromCol p)
                         in do
                             rowOper op1 mat
                             rowOper op2 mat
                             g <- foldM (\acc c -> fmap (fromIntegral . gcd acc . fromIntegral) $ readMatrix mat r c) 0 [p .. n-1]
                             if g == 0
                               then return()
                               else mapM_ (\c -> modifyMatrix mat r c (`quot` g)) [p .. n-1]

gaussian' :: Int -> Int -> Int -> Int -> STMatrix s Double -> ST s ()
gaussian' m n i j mat = do
    iPivot' <- findPivotMax m n i j mat
    case iPivot' of
        Nothing     -> return ()
        Just (r, p) -> do
                          rowOper (SWAP i r (FromCol j)) mat
                          pv <- readMatrix mat i p
                          mapM_ (reduce pv p) [i+1 .. m-1]
                          gaussian' m n (i+1) (p+1) mat
  where
    reduce pv p r = do
                      rv <- readMatrix mat r p
                      if abs rv < eps
                        then return ()
                        else
                         let frac = -rv / pv
                             op = AXPY frac i r (FromCol p)
                         in do
                             rowOper op mat
                             mapM_ (\j' -> modifyMatrix mat r j' (\x -> if abs x < eps then 0 else x)) [p..n-1]

-- | Gaussian elimination perfomed in-place in the @'ST'@ monad.

gaussianFFST :: Int -> Int -> STMatrix s Z -> ST s ()
gaussianFFST m n = gaussianFF' m n 0 0

gaussianST :: Int -> Int -> STMatrix s Double -> ST s ()
gaussianST m n = gaussian' m n 0 0


-- | Gaussian elimination as pure function. Involves a copy of the input matrix.
--
-- @
-- &#x3BB; let mat = (3 >< 4) [1, 1, -2, 0, 0, 2, -6, -4, 3, 0, 3, 1]
-- &#x3BB; mat
-- (3><4)
--  [ 1.0, 1.0, -2.0,  0.0
--  , 0.0, 2.0, -6.0, -4.0
--  , 3.0, 0.0,  3.0,  1.0 ]
-- &#x3BB; gaussian mat
-- (3><4)
--  [ 3.0, 0.0,  3.0,                1.0
--  , 0.0, 2.0, -6.0,               -4.0
--  , 0.0, 0.0,  0.0, 1.6666666666666667 ]
-- @
--

gaussianFF :: Matrix Z -> Matrix Z
gaussianFF mat = runST $ do
    matST <- thawMatrix mat
    gaussianFFST m n matST
    freezeMatrix matST
  where
    m = rows mat
    n = cols mat

gaussian :: Matrix Double -> Matrix Double
gaussian mat = runST $ do
    matST <- thawMatrix mat
    gaussianST m n matST
    freezeMatrix matST
  where
    m = rows mat
    n = cols mat

-- | Returns the indices of a maximal linearly independent subset of the columns
--   in the matrix.
--
-- @
-- &#x3BB; let mat = (3 >< 4) [1, 1, -2, 0, 0, 2, -6, -4, 3, 0, 3, 1]
-- &#x3BB; mat
-- (3><4)
--  [ 1.0, 1.0, -2.0,  0.0
--  , 0.0, 2.0, -6.0, -4.0
--  , 3.0, 0.0,  3.0,  1.0 ]
-- &#x3BB; independentColumns mat
-- [0,1,3]
-- @
--

independentColumnsFF :: Matrix Z -> [Int]
independentColumnsFF mat = pivotsUFF mat'
    where
        mat' = gaussianFF mat

independentColumns :: Matrix Double -> [Int]
independentColumns mat = pivotsU mat'
    where
        mat' = gaussian mat

-- | Returns a sub matrix containing a maximal linearly independent subset of
--   the columns in the matrix.
--
-- @
-- &#x3BB; let mat = (3 >< 4) [1, 1, -2, 0, 0, 2, -6, -4, 3, 0, 3, 1]
-- &#x3BB; mat
-- (3><4)
--  [ 1.0, 1.0, -2.0,  0.0
--  , 0.0, 2.0, -6.0, -4.0
--  , 3.0, 0.0,  3.0,  1.0 ]
-- &#x3BB; independentColumnsMat mat
-- (3><3)
--  [ 1.0, 1.0,  0.0
--  , 0.0, 2.0, -4.0
--  , 3.0, 0.0,  1.0 ]
-- @
--

independentColumnsMatFF :: Matrix Z -> Matrix Z
independentColumnsMatFF mat =
  case independentColumnsFF mat of
    [] -> (rows mat >< 1) $ repeat 0
    cs -> mat ¿ cs

independentColumnsMat :: Matrix Double -> Matrix Double
independentColumnsMat mat =
  case independentColumns mat of
    [] -> (rows mat >< 1) $ repeat 0
    cs -> mat ¿ cs
