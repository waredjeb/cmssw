import sys
import os
import torch

# prepare the datadir
if len(sys.argv) >= 2:
    datadir = sys.argv[1]
else:
    thisdir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.join(os.path.dirname(thisdir), "bin", "data")

os.makedirs(datadir, exist_ok=True)

class MultiTensorModel(torch.nn.Module):
    def __init__(self):
        super(MultiTensorModel, self).__init__()

        
    def forward(self, position, rotation, cluster, count):
        x = position * rotation[:, 0].unsqueeze(1) + position * rotation[:, 1].unsqueeze(1)
        x = x[:, :2]
        count = torch.cat([count, count], axis=1)
        x += count
        y = cluster[:, 0, :].unsqueeze(2) @ cluster[:, 0, :].unsqueeze(1)

        # TODO: Make multi output
        return x


module = MultiTensorModel()
position = torch.tensor([[0., 1., 2.], [1., 2., 3.], [2., 3., 4.], [3., 4., 5.]])
rotation = torch.tensor([[0., 0.], [1., 2.], [2., 4.], [3., 6.]])
cluster = torch.tensor([[[0., 1., 2.]], [[1., 2., 3.]], [[2., 3., 4.]], [[3., 4., 5.]], [[4., 5., 6.]], [[5., 6., 7.]]])
count = torch.tensor([[6], [6], [6], [6]])


print(module(position, rotation, cluster, count))
tm = torch.jit.trace(module.eval(), (position, rotation, cluster, count))
tm.save(f"{datadir}/multi_input_output_model.pt")

print("multi_input_output_model.pt created successfully!")