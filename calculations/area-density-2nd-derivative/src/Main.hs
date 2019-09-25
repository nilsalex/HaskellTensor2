{-#LANGUAGE DataKinds #-}
{-#LANGUAGE RankNTypes #-}

import Math.Tensor
import Math.Tensor.LorentzGenerator
import Math.Tensor.Examples.Gravity
import Math.Tensor.Examples.Gravity.DiffeoSymEqns
import Math.Tensor.Examples.Gravity.Schwarzschild

import Math.Tensor.Internal.LinearAlgebra

import Numeric.LinearAlgebra (rank)
import Numeric.LinearAlgebra.Data (toLists, fromLists, size, cmap, (===), Matrix, Z)
import Data.List (nub, nubBy, findIndex)

import Data.Maybe (mapMaybe)

import Data.Int (Int64)
import Data.Ratio

import Control.Parallel.Strategies (parList, runEvalIO, rdeepseq)

import qualified Data.IntMap.Strict as I

main = do

  let ans0 = fromListT6 [((Empty, Empty, Empty, Empty, Empty, Empty), AnsVar $ I.singleton 1 (SField 1))] :: ATens 0 0 0 0 0 0 AnsVarR
  let (eta4,eps4,ans4) = mkAnsatzTensorFastAbs 4 symList4 areaList4 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 1 0 0 0 0 0 AnsVarR)
  let (eta6,eps6,ans6) = mkAnsatzTensorFastAbs 6 symList6 areaList6 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 1 0 1 0 0 0 AnsVarR)
  let (eta8,eps8,ans8) = mkAnsatzTensorFastAbs 8 symList8 areaList8 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 2 0 0 0 0 0 AnsVarR)
  let (eta10_1,eps10_1,ans10_1) = mkAnsatzTensorFastAbs 10 symList10_1 areaList10_1 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 2 0 0 0 2 0 AnsVarR)
  let (eta10_2,eps10_2,ans10_2) = mkAnsatzTensorFastAbs 10 symList10_2 areaList10_2 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 2 0 1 0 0 0 AnsVarR)
  let (eta12,eps12,ans12) = mkAnsatzTensorFastAbs 12 symList12 areaList12 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 3 0 0 0 0 0 AnsVarR)
  let (eta14_1,eps14_1,ans14_1) = mkAnsatzTensorFastAbs 14 symList14_1 areaList14_1 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 3 0 0 0 2 0 AnsVarR)
  let (eta14_2,eps14_2,ans14_2) = mkAnsatzTensorFastAbs 14 symList14_2 areaList14_2 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 3 0 1 0 0 0 AnsVarR)

  let r0    = tensorRank6' ans0
  let r4    = tensorRank6' ans4
  let r6    = tensorRank6' ans6
  let r8    = tensorRank6' ans8
  let r10_1 = tensorRank6' ans10_1
  let r10_2 = tensorRank6' ans10_2
  let r12   = tensorRank6' ans12
  let r14_1 = tensorRank6' ans14_1
  let r14_2 = tensorRank6' ans14_2

  let ans0'    = ans0
  let ans4'    = shiftLabels6 r0 ans4
  let ans6'    = shiftLabels6 (r0 + r4) ans6
  let ans8'    = shiftLabels6 (r0 + r4 + r6) ans8
  let ans10_1' = shiftLabels6 (r0 + r4 + r6 + r8) ans10_1
  let ans10_2' = shiftLabels6 (r0 + r4 + r6 + r8 + r10_1) ans10_2
  let ans12'   = shiftLabels6 (r0 + r4 + r6 + r8 + r10_1 + r10_2) ans12
  let ans14_1' = shiftLabels6 (r0 + r4 + r6 + r8 + r10_1 + r10_2 + r12) ans14_1
  let ans14_2' = shiftLabels6 (r0 + r4 + r6 + r8 + r10_1 + r10_2 + r12 + r14_1) ans14_2

  let ansaetze = ans0' &.&>
                 ans4' &.&>
                 ans6' &.&>
                 ans8' &.&>
                 ans10_1' &.&>
                 ans10_2' &.&>
                 ans12' &.&>
                 ans14_1' &.&>
                 (singletonTList6 ans14_2')

  let e1 = eqn1 ans0' ans4'
  let e2 = eqn3 ans6'
  let e3 = eqn1A ans4' ans8'
  let e4 = eqn1AI ans6' ans10_2'
  let e5 = eqn2Aa ans6' ans10_1'
  let e6 = eqn3A ans6' ans10_2'
  let e7 = eqn1AB ans8' ans12'
  let e8 = eqn1ABI ans10_2' ans14_2'
  let e9 = eqn1AaBb ans10_1' ans14_1'
  let e10 = eqn2ABb ans10_1' ans10_2' ans14_1'
  let e11 = eqn3AB ans10_2' ans14_2'

  let system = e1 &.&>
               e2 &.&>
               e3 &.&>
               e4 &.&>
               e5 &.&>
               e6 &.&>
               e7 &.&>
               e8 &.&>
               e9 &.&>
               e10 &.&>
               (singletonTList6 e11)

  let mat' = map (map round) $ toLists $ toMatrixT6 system :: [[Int64]]
  let mat'' = filter (\rs -> any (/=0) rs) $ mat'
  let mat''' = nubBy compRows mat''
  let mat = fromLists $ mat'''

  let mat2 = rref mat

  let ps = pivotsUFF mat2

  let matD = cmap fromIntegral mat :: Matrix Double
  let mat2D = cmap fromIntegral mat2 :: Matrix Double

  print $ rank matD
  print $ rank mat2D
  print $ rank $ matD === mat2D
  print $ length ps

  --print $ (toLists mat2) !! 187

  --mapM_ print $ toLists mat2
  --mapM_ print mat'''
  let sol = fromRref mat2

  let ans0'' = applySolutionToTensor sol ans0'
  let ans4'' = applySolutionToTensor sol ans4'
  let ans6'' = applySolutionToTensor sol ans6'
  let ans8'' = applySolutionToTensor sol ans8'
  let ans10_1'' = applySolutionToTensor sol ans10_1'
  let ans10_2'' = applySolutionToTensor sol ans10_2'
  let ans12'' = applySolutionToTensor sol ans12'
  let ans14_1'' = applySolutionToTensor sol ans14_1'
  let ans14_2'' = applySolutionToTensor sol ans14_2'

  let solution = ans0'' &.&>
                 ans4'' &.&>
                 ans6'' &.&>
                 ans8'' &.&>
                 ans10_1'' &.&>
                 ans10_2'' &.&>
                 ans12'' &.&>
                 ans14_1'' &.&>
                 (singletonTList6 ans14_2'')

  print $ tensorRank6 solution

  let e1' = removeZeros6 $ eqn1 ans0'' ans4''
  let e2' = removeZeros6 $ eqn3 ans6''
  let e3' = removeZeros6 $ eqn1A ans4'' ans8''
  let e4' = removeZeros6 $ eqn1AI ans6'' ans10_2''
  let e5' = removeZeros6 $ eqn2Aa ans6'' ans10_1''
  let e6' = removeZeros6 $ eqn3A ans6'' ans10_2''
  let e7' = removeZeros6 $ eqn1AB ans8'' ans12''
  let e8' = removeZeros6 $ eqn1ABI ans10_2'' ans14_2''
  let e9' = removeZeros6 $ eqn1AaBb ans10_1'' ans14_1''
  let e10' = removeZeros6 $ eqn2ABb ans10_1'' ans10_2'' ans14_1''
  let e11' = removeZeros6 $ eqn3AB ans10_2'' ans14_2''

  let system' = e1' &.&>
                e2' &.&>
                e3' &.&>
                e4' &.&>
                e5' &.&>
                e6' &.&>
                e7' &.&>
                e8' &.&>
                e9' &.&>
                e10' &.&>
                (singletonTList6 e11')

  print e1'
  print e2'
  print e3'
  print e4'
  print e5'
  print e6'
  print e7'
  print e8'
  print e9'
  print e10'
  print e11'

  return ()

type Solution = I.IntMap AnsVarR

fromRref :: Matrix Z -> Solution
fromRref ref = I.fromList assocs
    where
        rows   = toLists ref
        assocs = mapMaybe fromRow rows

fromRow :: Integral a => [a] -> Maybe (Int, AnsVarR)
fromRow xs = case assocs of
               []             -> Nothing
               [(i, _)]       -> Just (i, AnsVar I.empty)
               (i, v):assocs' -> let assocs'' = map (\(i, v') -> (i, SField $ - (fromIntegral v') % (fromIntegral v))) assocs'
                                 in Just (i, AnsVar $ I.fromList assocs'')
    where
        assocs = filter ((/=0) . snd) $ zip [1..] xs

applySolution :: Solution -> AnsVarR -> AnsVarR
applySolution s (AnsVar m) = AnsVar $
                             I.foldlWithKey' (\m' i sub -> let AnsVar m'' = AnsVar m' `addS` sub
                                                           in I.delete i m'') m s'
    where
        s' = I.intersectionWith (\row v -> v `prod` row) s m

applySolutionToTensor :: forall n1 n2 n3 n4 n5 n6.Solution -> ATens n1 n2 n3 n4 n5 n6 AnsVarR -> ATens n1 n2 n3 n4 n5 n6 AnsVarR
applySolutionToTensor sol t = removeZeros6 $ mapTo6 (applySolution sol) t

compRows :: (Integral a, Eq a) => [a] -> [a] -> Bool
compRows xs ys
        | ix /= iy  = False
        | x  == y   = xs == ys
        | otherwise = xs' == ys'
    where
        Just ix = findIndex (/= 0) xs
        Just iy = findIndex (/= 0) ys
        x = xs !! ix
        y = ys !! iy
        f = x%y
        xs' = map fromIntegral xs
        ys' = map ((*f) . fromIntegral) ys
