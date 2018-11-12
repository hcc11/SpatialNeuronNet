 function Col = HF_whiten(Col,whitePart)
   Col = Col + whitePart*([1,1,1]-Col);