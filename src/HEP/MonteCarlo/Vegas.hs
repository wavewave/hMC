{-# LANGUAGE BangPatterns, PackageImports #-}

module HEP.MonteCarlo.Vegas where

import System.Random.Mersenne

import Data.Array

import Control.Monad.State
import Control.Monad.IO.Class




data SamplingParameter2D = SP2D { 
  xStart  :: Double,  
  xEnd    :: Double, 
  xNumBin :: Int, 
  yStart  :: Double,
  yEnd    :: Double, 
  yNumBin :: Int, 
  numSamplePerBin :: Int
}

locateBin :: Double -> Double -> Int -> Double -> Int 
locateBin start end numbin x = 
  let step = (end - start) / fromIntegral numbin 
      nx = (x - start) / step
  in floor nx 

binRange start end numbin i = 
  let idouble = fromIntegral i
      step = (end - start) / (fromIntegral numbin)
  in  (start+idouble*step,start+(idouble+1)*step) 

 
data XorY = X | Y 

vegas2DSamplingM :: (MonadIO m) => 
                    ((Double,Double)-> m Double) -- ^ target function
                    -> SamplingParameter2D    -- ^ sampling parameter
                    -> MTGen                  -- ^ Random generator
                    -> m ((Array Int Double, Array Int Double), 
                          (Array Int Double, Array Int Double)) -- ^ Xpdf, Ypdf, 
                                                                --   Xcdf, Ycdf
vegas2DSamplingM f p gen = do 
  let num = numSamplePerBin p
  let sampleN (xi,xf) (yi,yf) n acc 
        | n <= 0 = return acc
        | n >  0 = do r1 <- liftIO $ random gen 
                      r2 <- liftIO $ random gen 
                      let x = r1*(xf-xi)+xi
                          y = r2*(yf-yi)+yi 
                      fxy <- f (x,y)    
                      sampleN (xi,xf) (yi,yf) (n-1) (acc + fxy)
  let sampleForBin xory i = do 
        let (rangex,rangey) = case xory of 
              X -> (binRange (xStart p) (xEnd p) (xNumBin p) i, (yStart p,yEnd p))
              Y -> ((xStart p,xEnd p), binRange (yStart p) (yEnd p) (yNumBin p) i) 
        r <- sampleN rangex rangey num 0.0 
--        putStrLn $ show r
        return $ r/(fromIntegral num)
        
  xsmpllst <- mapM (sampleForBin X) [0..xNumBin p - 1]
  ysmpllst <- mapM (sampleForBin Y) [0..yNumBin p - 1]
  
  let xsmparr = listArray (0,xNumBin p - 1) xsmpllst 
      ysmparr = listArray (0,yNumBin p - 1) ysmpllst
   
      (xpdf,xcdf) = cdfFromPdf xsmparr
      (ypdf,ycdf) = cdfFromPdf ysmparr      
      
  return ((xpdf,ypdf),(xcdf,ycdf))
  

vegasMCIntegrationM :: (MonadIO m) => 
                       ((Double,Double) -> m Double)    -- ^ function
                       -> SamplingParameter2D         -- ^ parameter
                       -> (Array Int Double, Array Int Double) -- ^ xpdf ypdf
                       -> (Array Int Double, Array Int Double) -- ^ xcdf ycdf
                       -> Int     -- ^ Number of monte carlo samples
                       -> MTGen
                       -> m Double 
vegasMCIntegrationM f p (xpdf,ypdf) (xcdf,ycdf) num gen = do  
  let worker !n !acc | n <= 0 = return acc  
                     | n >  0 = do 
                       x <- liftIO $ randomGenWithCDF xcdf (xStart p, xEnd p) gen
                       y <- liftIO $ randomGenWithCDF ycdf (yStart p, yEnd p) gen
                       let bx = locateBin (xStart p) (xEnd p) (xNumBin p) x
                           gxval = xpdf ! bx
                           by = locateBin (yStart p) (yEnd p) (yNumBin p) y
                           gyval = ypdf ! by
                       fxy <- f (x,y)
                       worker (n-1) $ acc + (fxy / gxval / gyval)
  r <- worker num 0.0
  return (r / fromIntegral (num*(xNumBin p)*(yNumBin p)))
                         
  
cdfFromPdf :: Array Int Double            -- ^ original unnormalized pdf
              -> (Array Int Double, Array Int Double)  -- ^ result normalized pdf and  cdf  
cdfFromPdf pdf = 
  let arraybound = bounds $ pdf 
      arrayboundPlus1 = (fst arraybound +1, snd arraybound +1)
      num = (+1) . snd $ arraybound 
      work n 
        | n == 0 = return ()
        | n > 0 = do (i,acc,lst) <- get 
                     let acc' = acc + pdf ! i
                     put (i+1, acc', acc':lst)  
                     work $ n-1
      (_,acc,revlst) = (flip execState) (0,0.0,[]) $ work num

      revdiv (x:xs) newlist d = revdiv xs ((x / d) : newlist) d
      revdiv [] newlist d = newlist

      lst = revdiv revlst [] acc
  in  (fmap (/acc) pdf, listArray arrayboundPlus1 lst)
                  
findBinFromCDF :: Array Int Double    -- ^ CDF      
                  -> Int              -- ^ arraybound
                  -> Double           -- ^ num
                  -> (Int,Int)  -- ^ (start,end)
                  -> Maybe Int              -- ^ result bin
findBinFromCDF cdf b num (start,end) 
  | num < 0.0 || num > 1.0 = Nothing 
  | start == end = Just start 
  | start == end -1 = Just start
  | otherwise    = let half = (start+end) `div` 2
                   in  if num > (cdf ! half) 
                       then findBinFromCDF cdf b num (half,end) 
                       else findBinFromCDF cdf b num (start,half)
                 
xformNumUsingCDF :: Array Int Double   -- ^ CDF 
                    -> (Double,Double) -- ^ Range
                    -> Double          -- ^ inputnum 
                    -> Double          -- ^ outputnum
xformNumUsingCDF cdf (start,end) num = do
  let (_,b) = bounds cdf
      bin2num i | i == 0 = 0.0 
                | otherwise = cdf ! i
  case findBinFromCDF cdf b num (0,b) of 
    Just bin -> let (binstart,binend) = binRange start end b bin 
                    ratio = (num - bin2num bin) / (bin2num (bin+1) - bin2num bin)
                in  (binend - binstart) * ratio + binstart
    Nothing  -> error "outofarray" -- putStrLn $ "outofarray"

randomGenWithCDF :: Array Int Double 
                    -> (Double,Double) 
                    -> MTGen
                    -> IO Double
randomGenWithCDF cdf range gen = do    
  r <- random gen 
  return $ xformNumUsingCDF cdf range r
  
                 