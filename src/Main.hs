{-# LANGUAGE DataKinds #-}

module Main (
 main
) where

import PerturbationTree2_3
import TensorTreeNumeric4_2
import BasicTensors4_2
import Data.List (intersperse, sortOn)
import Data.Ratio ((%), numerator, denominator)
import qualified Data.IntMap.Strict as I (fromList)

diffeo1 :: ATens 0 2 0 0 0 0 (AnsVar Rational) -> ATens 0 1 0 0 1 1 (AnsVar Rational)
diffeo1 ans8 = contrATens1 (0,1) $ ans8 &* flatInter

-- this equation differs from the scalar equations, block 2 is a density term!!
diffeo2 :: ATens 0 1 1 0 0 0 (AnsVar Rational) -> ATens 0 2 1 0 0 0 (AnsVar Rational) -> ATens 0 1 1 0 1 1 (AnsVar Rational)
diffeo2 ans6 ans10 = block1 &+ block2 &+ block3
    where
        block1 = contrATens1 (0,0) $ ans10 &* flatInter
        block2 = ans6 &* delta3
        block3 = contrATens1 (0,0) $ contrATens2 (0,0) $ ans6 &* interEqn3

diffeo3 :: ATens 0 1 1 0 0 0 (AnsVar Rational) -> ATens 0 2 0 0 2 0 (AnsVar Rational) -> ATens 0 1 0 0 3 1 (AnsVar Rational)
diffeo3 ans6 ans10 = block1 &+ block2
    where
        block1 = symATens5 (1,2) $ contrATens1 (0,1) $ ans10 &* flatInter
        block2 = contrATens1 (0,0) $ contrATens2 (0,0) $ 2 &. ans6 &* interEqn4

diffeo4 :: ATens 0 1 1 0 0 0 (AnsVar Rational) -> ATens 0 2 1 0 0 0 (AnsVar Rational) -> ATens 0 1 0 0 3 1 (AnsVar Rational)
diffeo4 ans6 ans10 = block1 &+ block2
    where
        block1 = contrATens2 (0,0) $ contrATens1 (0,1) $ contrATens1 (1,2) $ ans10 &* interEqn5 &* flatArea
        block2 = contrATens2 (0,0) $ contrATens1 (0,0) $ ans6 &* interEqn5

diffeo5 :: ATens 0 1 1 0 0 0 (AnsVar Rational) -> ATens 0 0 0 0 3 1 (AnsVar Rational)
diffeo5 ans6 = block1
    where
        block1 = contrATens2 (0,0) $ contrATens1 (0,0) $ contrATens1 (1,1) $ ans6 &* interEqn5 &* flatArea

diffeo6 :: ATens 0 0 0 0 0 0 (AnsVar Rational) -> ATens 0 0 0 0 1 1 (AnsVar Rational)
diffeo6 ans0 = delta3 &* ans0

alphabet :: String
alphabet = "abcdefghijklmnpqrstuvwxyz"

etaString :: Eta -> String
etaString (Eta a b) = "η_{" ++ [alphabet !! (a-1), alphabet !! (b-1)] ++ "}"

epsString :: Epsilon -> String
epsString (Epsilon a b c d) = "ε_{" ++ [alphabet !! (a-1), alphabet !! (b-1), alphabet !! (c-1), alphabet !! (d-1)] ++ "}"

etaList :: AnsatzForestEta -> [String]
etaList = map (\(f, x, s) -> show f ++ " e_" ++ show x ++ " " ++ s) .
          sortOn (\(_,x,_) -> x) .
          map (\(is, Var f x) -> (f, x, (concat $ intersperse " " $ map etaString is))) .
          flattenForest

epsList :: AnsatzForestEpsilon -> [String]
epsList = map (\(f, x, s) -> show f ++ " e_" ++ show x ++ " " ++ s) .
          sortOn (\(_,x,_) -> x) .
          map (\(e, is, Var f x) -> (f, x, (concat $ intersperse " " $ [epsString e] ++ map etaString is))) .
          flattenForestEpsilon


invEta2 :: ATens 0 0 2 0 0 0 Rational
invEta2 = contrATens3 (0,0) $ contrATens3 (1,1) $ contrATens3 (1,2) $ contrATens3 (3,3) $
          interI2Inv &* interI2 &* invEta &* invEta

main = do 
    --mass term ansätze
    
    let (ans8ABEta,ans8ABEps,ans8AB') = mkAnsatzTensorFast 8 filterList8 symList8 areaList8IndsEta areaList8IndsEps 

    --kinetic term ansätze 

    let (ans6AIEta,ans6AIEps,ans6AI'') = mkAnsatzTensorFast 6 filterList6 symList6 areaList6IndsEta areaList6IndsEps
    let (ans10ApBqEta,ans10ApBqEps,ans10ApBq'') = mkAnsatzTensorFast 10 filterList10_1 symList10_1 areaList10_1IndsEta areaList10_1IndsEps
    let (ans10ABIEta,ans10ABIEps,ans10ABI'') = mkAnsatzTensorFast 10 filterList10_2 symList10_2 areaList10_2IndsEta areaList10_2IndsEps

    let ans6AI'    = contrATens2 (1,0) $ invEta2 &* ans6AI''
    let ans10ApBq' = contrATens3 (1,0) $ contrATens3 (3,1) $ invEta &* invEta &* ans10ApBq''
    let ans10ABI'  = contrATens2 (1,0) $ invEta2 &* ans10ABI''

    let r6      = tensorRank' ans6AI'
    let r8      = tensorRank' ans8AB'
    let r10ApBq = tensorRank' ans10ApBq'
    let r10ABI  = tensorRank' ans10ABI'
    let r = r6 + r8 + r10ApBq + r10ABI

    let ans8AB    = ans8AB'                                      -- from 1 to 6
    let ans10ApBq = shiftLabels6 r8 ans10ApBq'                   -- from 7 to 21
    let ans10ABI  = shiftLabels6 (r8 + r10ApBq) ans10ABI'        -- from 22 to 37
    let ans6AI    = shiftLabels6 (r8 + r10ApBq + r10ABI) ans6AI' -- from 38 to 40

    let d1 = diffeo1 ans8AB
    let d2 = diffeo2 ans6AI ans10ABI
    let d3 = diffeo3 ans6AI ans10ApBq
    let d4 = diffeo4 ans6AI ans10ABI
    let d5 = diffeo5 ans6AI

    let system = d5 &> d4 &> d3 &> d2 &> (singletonTList d1)
--    let system = singletonTList $ diffeo1 ans8AB'

--    putStrLn $ "DOFs      : " ++ (show r)
--    putStrLn $ "my eqns   : " ++ (show $ tensorRank system)

--    sequence_ $ map putStrLn $ etaList ans10ApBqEta
--    sequence_ $ map putStrLn $ epsList ans10ApBqEps
    putStrLn $ unlines $ map (\((i, j), v) -> if denominator v /= 1
                                              then undefined
                                              else "(" ++ show i ++ ", " ++ show j ++ ") = " ++ show (numerator v) ++ ",")
                       $ toMatList6 system
