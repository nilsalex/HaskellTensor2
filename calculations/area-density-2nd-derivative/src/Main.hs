{-# LANGUAGE DataKinds #-}

import Math.Tensor
import Math.Tensor.LorentzGenerator
import Math.Tensor.Examples.Gravity
import Math.Tensor.Examples.Gravity.DiffeoSymEqns
import Math.Tensor.Examples.Gravity.Schwarzschild

import Math.Tensor.Internal.LinearAlgebra

import Numeric.LinearAlgebra (rank)
import Numeric.LinearAlgebra.Data (toLists, fromLists, size, cmap, (===), Matrix, Z)
import Data.List (nub, nubBy, findIndex, sort)

import Data.Maybe (mapMaybe)

import Data.Int (Int64)
import Data.Ratio

import Control.Parallel.Strategies (parList, runEvalIO, rdeepseq)

import qualified Data.IntMap.Strict as I

main = do

  let ans0 = fromListT6 [((Empty, Empty, Empty, Empty, Empty, Empty), AnsVar $ I.singleton 1 (SField 1))] :: ATens 0 0 0 0 0 0 AnsVarR
--  let (eta4,eps4,ans4) = mkAnsatzTensorFastAbs 4 symList4 areaList4 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 0 1 0 0 0 0 AnsVarR)
  let ans4 = ZeroTensor :: ATens 0 1 0 0 0 0 AnsVarR
  let (eta6,eps6,_ans6) = mkAnsatzTensorFastAbs 6 symList6 areaList6 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 0 1 0 1 0 0 AnsVarR)
  let (eta8,eps8,_ans8) = mkAnsatzTensorFastAbs 8 symList8 areaList8 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 0 2 0 0 0 0 AnsVarR)
  let (eta10_1,eps10_1,_ans10_1) = mkAnsatzTensorFastAbs 10 symList10_1 areaList10_1 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 0 2 0 0 0 2 AnsVarR)
  let (eta10_2,eps10_2,_ans10_2) = mkAnsatzTensorFastAbs 10 symList10_2 areaList10_2 :: (AnsatzForestEta, AnsatzForestEpsilon, ATens 0 2 0 1 0 0 AnsVarR)

  let met = fromListT6 $ map (\i -> let s = case i of
                                              0 -> 1
                                              1 -> -1
                                              2 -> -1
                                              3 -> -1
                                              4 -> 1
                                              5 -> 1
                                              6 -> 1
                                              7 -> 1
                                              8 -> 1
                                              9 -> 1
                                              _  -> error "Ind9 greater than 9"
                                    in ((Empty, Empty, (Ind9 i) `Append` ((Ind9 i) `Append` Empty), Empty, Empty, Empty), SField s))
                         [0..9] :: ATens 0 0 2 0 0 0 (SField Rational)

  let ans6 = contrATens2 (0,0) $ met &* _ans6
  let ans8 = _ans8
  let ans10_1 = contrATens3 (0,0) $ contrATens3 (2,1) $ invEtaA &* invEtaA &* _ans10_1
  let ans10_2 = contrATens2 (0,0) $ met &* _ans10_2

  let r0    = tensorRank6' ans0
  let r4    = tensorRank6' ans4
  let r6    = tensorRank6' ans6
  let r8    = tensorRank6' ans8
  let r10_1 = tensorRank6' ans10_1
  let r10_2 = tensorRank6' ans10_2

  let ans0'    = ans0
  let ans4'    = shiftLabels6 r0 ans4
  let ans6'    = shiftLabels6 (r0 + r4) ans6
  let ans8'    = shiftLabels6 (r0 + r4 + r6) ans8
  let ans10_1' = shiftLabels6 (r0 + r4 + r6 + r8) ans10_1
  let ans10_2' = shiftLabels6 (r0 + r4 + r6 + r8 + r10_1) ans10_2

  let ansaetze = ans0' &.&>
                 ans4' &.&>
                 ans6' &.&>
                 ans8' &.&>
                 ans10_1' &.&>
                 (singletonTList6 ans10_2')

  putStrLn $ "dimension of ansatz space   : " ++ show (tensorRank6 ansaetze)

  let two = SField (2 :: Rational)

  let e1 = eqn1 ans0' ans4'
  let e2 = eqn3 ans6'
  let e3 = eqn1A ans4' (two &. ans8')
  let e4 = eqn1AI ans6' ans10_2'
  let e5 = eqn2Aa ans6' (two &. ans10_1')
  let e6 = eqn3A ans6' ans10_2'

  let system = e1 &.&>
               e2 &.&>
               e3 &.&>
               e4 &.&>
               e5 &.&>
               (singletonTList6 e6)

  putStrLn $ "rank of system              : " ++ show (tensorRank6 system)

  let solution = redefineVarsSystem6 $ solveSystem6 system ansaetze

  case solution of
    (t1 `AppendTList6` (
     t2 `AppendTList6` (
     t3 `AppendTList6` (
     t4 `AppendTList6` (
     t5 `AppendTList6` (
     t6 `AppendTList6` (
     EmptyTList6)))))))
      -> do
          let Just ans0'' = tryAsATens ans0 t1
          let Just ans4'' = tryAsATens ans4 t2
          let Just ans6'' = tryAsATens ans6 t3
          let Just ans8'' = tryAsATens ans8 t4
          let Just ans10_1'' = tryAsATens ans10_1 t5
          let Just ans10_2'' = tryAsATens ans10_2 t6

          let e1' = removeZeros6 $ eqn1 ans0'' ans4''
          let e2' = removeZeros6 $ eqn3 ans6''
          let e3' = removeZeros6 $ eqn1A ans4'' (two &. ans8'')
          let e4' = removeZeros6 $ eqn1AI ans6'' ans10_2''
          let e5' = removeZeros6 $ eqn2Aa ans6'' (two &. ans10_1'')
          let e6' = removeZeros6 $ eqn3A ans6'' ans10_2''

          putStrLn $ "dimension of solution space : " ++ show (tensorRank6 solution)
          putStrLn ""
          putStrLn "Equations on solution space :"

          print e1'
          print e2'
          print e3'
          print e4'
          print e5'
          print e6'

          let kin' = ans10_2'' &- (contrATens3 (0,0) $ contrATens3 (1,1) $ interI2 &* ans10_1'')
          let kin = kin' &+ (tensorTrans2 (0,1) kin')

          putStrLn ""
          putStrLn $ "number of parameters in kinetic part of EOM : " ++ show (tensorRank6' kin)
          putStrLn $ "number of parameters in mass part of EOM    : " ++ show (tensorRank6' ans8'')
 
          let n1 = removeZeros6 $ (delta3A &* ans4'') &+ (contrATens1 (0,1) $ interArea &* ans4'') &+ (contrATens1 (0,0) $ contrATens1 (1,1) $ two &. flatArea &* interArea &* ans8'')

          let n2_1 = contrATens1 (0,0) $ contrATens1 (1,1) $ contrATens2 (0,0) $ flatArea &* interArea &* interJ2 &* ans10_2''
          let n2_2 = contrATens1 (0,0) $ contrATens1 (1,2) $ contrATens2 (0,0) $ flatArea &* interArea &* interJ2 &* ans10_2''
          let n2_3 = contrATens1 (0,0) $ contrATens1 (1,1) $ flatArea &* interArea &* ans10_1''

          let n2 = removeZeros6 $ cyclicSymATens5 [0,1,2] $ n2_1 &+ n2_2 &- (two &. n2_3)

          putStrLn ""
          putStrLn "Noether identity 1 :"
          print n1

          putStrLn ""
          putStrLn "Noether identity 2 :"
          print n2

          print $ toListShow6 ans0
          print $ toListShow6 ans4
          print $ toListShow6 ans0''
          print $ toListShow6 ans4''

  return ()
