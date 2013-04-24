{-# LANGUAGE BangPatterns #-}

-----------------------------------------------------------------------------
-- |
-- Module      : Numeric.MonteCarlo.Plain
-- Copyright   : (c) 2011, 2013 Ian-Woo Kim
--
-- License     : GPL-3
-- Maintainer  : Ian-Woo Kim <ianwookim@gmail.com>
-- Stability   : experimental
-- Portability : GHC
--
-- Plain non-adaptive monte carlo integration algorithm  
--
-----------------------------------------------------------------------------

module Numeric.MonteCarlo.Plain where

import System.Random.Mersenne

-- | 1D integration
plain1DMC :: (Double -> IO Double) 
             -> Int 
             -> IO Double 
plain1DMC integrand num  = do  
  g <- getStdGen 
  let go :: Double -> Int -> IO Double 
      go !accr !accn   
        | accn <= 0 = 
          return accr
        | otherwise = do 
          x1 <- random g :: IO Double 
          itr <- integrand x1 
          if isNaN itr 
            then putStrLn $ show (x1,itr)
            else return ()
          go (accr + itr) (accn-1) 
      
  r <- go 0.0 num 
  return $ r/fromIntegral num 


-- | 2D integration
plain2DMC :: ((Double,Double) -> IO Double) 
             -> Int 
             -> IO Double 
plain2DMC integrand num  = do  
  g <- getStdGen 
  let go :: Double -> Int -> IO Double 
      go !accr !accn   
        | accn <= 0 = 
          return accr
        | otherwise = do 
          x1 <- random g :: IO Double 
          x2 <- random g :: IO Double 
          itr <- integrand (x1,x2) 
          if isNaN itr 
            then putStrLn $ show (x1,x2,itr)
            else return ()
          go (accr + itr) (accn-1) 
      
  r <- go 0.0 num 
  return $ r/fromIntegral num 

-- | 3D integration
plain3DMC :: ((Double,Double,Double) -> IO Double) 
             -> Int 
             -> IO Double 
plain3DMC integrand num  = do  
  g <- getStdGen 
  let go :: Double -> Int -> IO Double 
      go !accr !accn   
        | accn <= 0 = 
          return accr
        | otherwise = do 
          x1 <- random g :: IO Double 
          x2 <- random g :: IO Double 
          x3 <- random g :: IO Double 
          itr <- integrand (x1,x2,x3) 
          if isNaN itr 
            then putStrLn $ show (x1,x2,x3,itr)
            else return ()
          go (accr + itr) (accn-1) 
      
  r <- go 0.0 num 
  return $ r/fromIntegral num 


