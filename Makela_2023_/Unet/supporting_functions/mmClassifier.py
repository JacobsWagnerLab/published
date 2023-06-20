import torch
from torch.autograd import Variable
import torch.nn.functional as F
from torch.utils.data import DataLoader
from tqdm import tqdm
import numpy as np
from skimage.measure import label
from skimage import io
import matplotlib.pyplot as plt

class mmClassifier:
    
    """ trainer for microscopy experiment UNet. 
    The classifier will run 1 epoch, then test the accuracy on test data."""
    
    def __init__(self,net,optimizer,scheduler,threshold,num_epochs,save_name):
        """
        The classifier for mm experiment
            net (nn.Module): The neural net module containing the definition of your model
            num_epochs (int): The maximum number of epochs on which the model will train
        """
        self.net = net
        self.epoch_counter = 0
        self.use_cuda = torch.cuda.is_available()
        self.num_epochs = num_epochs
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.threshold = threshold
        self.save_name = save_name

            
    def _criterion(self, pred, target):
        #defime loss as BCE + DC
        
        #release from logits
        pred = torch.sigmoid(pred)

        num = target.size(0)
        m1 = pred.view(num,-1)
        m2 = target.view(num,-1)
        
        # binary cross entropy 
        BCE = F.binary_cross_entropy(m1,m2)
        
        #dice loss
        intersection = (m1 * m2)
        DC = 2. * (intersection.sum(1)) / (m1.sum(1) + m2.sum(1))
        DC = 1 - DC.sum() / num

        return DC + 0.5*BCE
        
        #return F.binary_cross_entropy(segmentation, target)
    
    def _trainEpoch(self, train_loader):
        #run training on 1 batch
        print('training epoch {}; training losses:'.format(self.epoch_counter+1))
        running_losses = []
        # with tqdm(total = len(train_loader)) as pbar:
        for i_batch, sample_batched in enumerate(train_loader):
            im, target = sample_batched['im'], sample_batched['mask']
            im = im.cuda() 
            self.optimizer.zero_grad()
            loss=0

            output = self.net(im)

            # if deep supervision for nested Unet
            if len(output) > 1:
                for res in output:
                    loss += self._criterion(res, target.cuda())
                loss /= len(output)
            else: 
                loss = self._criterion(output, target.cuda())
            
            running_losses.append(loss.item())
            loss.backward()
            self.optimizer.step()
            # pbar.update(1)

        print("epoch {} loss: {} /n".format(self.epoch_counter+1, np.mean(running_losses)))        
        # pbar.close()
                    
    def _validateEpoch(self, validation_loader):
        
        #run network on validation data after each epoch.
        #TODO - select random slices form the dataset, save them in a folder, measure accuracies
        accuracy = []        
        for i_batch, sample_batched in enumerate(validation_loader):
            im,_ = sample_batched['im'], sample_batched['mask']
            im = im.cuda() 
            res = self.net(im)
            
            #if deep supervision for nested Unet
            if len(res) > 1:
                res = res[-1]

            res = torch.sigmoid(res)
            res = res.to("cpu").detach().numpy().squeeze(0).squeeze(0)
            #res = res > self.threshold
            #res = res.astype('uint16')
            
            #target =  target.to("cpu").detach().numpy().squeeze(0).squeeze(0).astype('uint16')
            
            file_name_seg = 'epoch_' + str(self.epoch_counter+1) + '_' + str(i_batch) + '_seg.tiff'
            #file_name_target = 'epoch_' + str(self.epoch_counter) + '_' + str(i_batch) + '_target.tiff'
            io.imsave('/hdd/RecPAIR/PraNetTraining/' + file_name_seg, res, compress=6)
            #io.imsave('/hdd/RecPAIR/PraNetTraining/' + file_name_target, target, compress=6)

            if i_batch > 0:
                torch.cuda.empty_cache()
                break

        #accuracy.append(np.abs(np.max(res_labels) - np.max(target_labels)) /  np.max(target_labels))
        
        
    def _runEpoch(self, train_loader: DataLoader, validation_loader: DataLoader, callbacks=None):
        # run a single epoch and validate the output. Possibly save it too
        
        self.net.train()

        # Run a train pass on the current epoch
        self._trainEpoch(train_loader)

        # switch to evaluate mode
        self.net.eval()

        # Run the validation pass
        # self._validateEpoch(validation_loader)
        
    def train(self, train_loader: DataLoader, validation_loader: DataLoader):
        
       # if self.use_cuda:
       #     self.net.cuda()
        
        #call training of the network
        for epoch in range(self.num_epochs):
            self._runEpoch(train_loader, validation_loader)
            self.scheduler.step()
            self.epoch_counter = self.epoch_counter + 1
        
        torch.save({
                 'epoch': self.epoch_counter,
                 'model_state_dict' : self.net.state_dict(),
                 'optimizer_state_dict' : self.optimizer.state_dict()
             },  self.save_name)
             