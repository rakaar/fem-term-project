function [nodeCoordinates,elementNodes] = rectangularMesh(Lx,Ly, numberElementsX,numberElementsY)  

xx = 0:Lx/numberElementsX:Lx;
 yy = 0:Ly/numberElementsY:Ly;
 [XX YY] = meshgrid(yy,xx);
 nodeCoordinates = [YY(:),XX(:)];
 elementNodes = zeros(numberElementsX*numberElementsY,4);  j = 1;
 i =1;
 i1 =0;
 counter = 0;
 for j = 1:numberElementsY
 for i =1: numberElementsX
 counter = counter +1;
 if i ==1 && j==1
 i1 =1;
else
 i1 = i1 +1;
end

 i2 = i1 + 1;
 i4 = i1 + numberElementsX + 1;
 i3 = i2 + numberElementsX + 1;
 elementNodes(counter,:) = [i1 i2 i3 i4]; 
 end
 i1 = i1+1;i2 = i2+1;
 end
end