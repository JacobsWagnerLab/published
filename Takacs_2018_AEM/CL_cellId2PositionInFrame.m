function positionInFrame = CL_cellId2PositionInFrame(cellId, frame, CL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function positionInFrame = CL_cellId2PositionInFrame(cellId, frame, CL)
%microbeTracker.v0.2.8
%@author:  Ahmad J Paintdakhi
%@date:    March 21, 2013
%@copyright 2012-2013 Yale University
%==========================================================================
%**********output********:
%positionInFrame:   position of given cellId in a given frame.
%**********Input********:
%CL:        A structure containing two fields meshData and cellId
%frame:     frame number
%cellId:    cell id in a given frame.
%==========================================================================
%The function returns the position of cellId in a given frame.
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
positionInFrame = [];
if CL_isNonEmptyFrame(frame, CL)
   positionInFrame = find(CL.cellId{frame}==cellId);
end
end