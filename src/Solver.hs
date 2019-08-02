module Solver where

import qualified Data.Map.Strict as M

type Matrix a = [[a]]
type Row a = [a]

fromSparse :: Num a => Int -> Int -> M.Map (Int, Int) a -> [[a]]
fromSparse rows cols m =
    (map . map) (\(r,c) -> case M.lookup (r, c) m
                           of Nothing -> 0
                              Just v  -> v)
    numbered
  where
    numbered = map (\r -> map (\c -> (r, c)) [0..cols-1]) [0..rows-1]

pivotInRow :: (Num a, Eq a) => [a] -> Maybe Int
pivotInRow [] = Nothing
pivotInRow (x:xs) = if x == 0
                    then (+) <$> (Just 1) <*> pivotInRow xs
                    else Just 0

incrementFirstMaybe :: Maybe (Int, Int) -> Maybe (Int, Int)
incrementFirstMaybe Nothing = Nothing
incrementFirstMaybe (Just (r, c)) = Just (r+1, c)

pivotInMat :: (Num a, Eq a) => [[a]] -> Maybe (Int, Int)
pivotInMat [] = Nothing
pivotInMat (r:rs) = case pivotInRow r
                    of Nothing -> incrementFirstMaybe (pivotInMat rs)
                       Just pr -> Just (0, pr)

swap :: Int -> Int -> [a] -> [a]
swap _ _ [] = error "cannot swap in empty list"
swap _ _ (x:[]) = error "cannot swap in singleton list"
swap a b xs@(x:xs')
    | a < 0 || b < 0 = error "wrong parameters in swap"
    | a == b = xs
    | b < a = swap b a xs
    | a > 0 = x : (swap (a-1) (b-1) xs')
    | a == 0 = let (l1, l2) = splitAt b xs'
               in (last l1) : (init l1 ++ [x] ++ l2)

normalize :: Fractional a => Int -> [a] -> [a]
normalize p xs = map (/pv) xs
    where
        pv = xs !! p
