{-# LANGUAGE BangPatterns #-}

module HEP.MonteCarlo.Plain where

import System.Random.Mersenne

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


