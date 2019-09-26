{-# LANGUAGE DataKinds #-}

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

  putStrLn $ "dimension of ansatz space   : " ++ show (tensorRank6 ansaetze)

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

  putStrLn $ "rank of system              : " ++ show (tensorRank6 system)

  let solution = solveSystem6 system ansaetze

  case solution of
    (t1 `AppendTList6` (
     t2 `AppendTList6` (
     t3 `AppendTList6` (
     t4 `AppendTList6` (
     t5 `AppendTList6` (
     t6 `AppendTList6` (
     t7 `AppendTList6` (
     t8 `AppendTList6` (
     t9 `AppendTList6` (
     EmptyTList6))))))))))
      -> do
          let Just ans0'' = tryAsATens ans0 t1
          let Just ans4'' = tryAsATens ans4 t2
          let Just ans6'' = tryAsATens ans6 t3
          let Just ans8'' = tryAsATens ans8 t4
          let Just ans10_1'' = tryAsATens ans10_1 t5
          let Just ans10_2'' = tryAsATens ans10_2 t6
          let Just ans12'' = tryAsATens ans12 t7
          let Just ans14_1'' = tryAsATens ans14_1 t8
          let Just ans14_2'' = tryAsATens ans14_2 t9

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

          putStrLn $ "dimension of solution space : " ++ show (tensorRank6 solution)
          putStrLn ""
          putStrLn "Equations on solution space :"

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
