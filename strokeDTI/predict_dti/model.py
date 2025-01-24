import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.nn import Linear, Parameter
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree
from torch_geometric.data import Data
from torch_scatter import scatter_mean
from torch_geometric.nn import (
    GraphConv,
    global_add_pool,
    TransformerConv,
    GATv2Conv,
    GENConv,
    GeneralConv,
    GINEConv,
    ResGatedGraphConv,
)
from torch_scatter import scatter_mean

from strokeDTI.predict_dti.params import DEVICE


def get_model_from_name(model_name):
    if model_name == "transformer_cnn":
        model = Transformer()

    elif model_name == "gatv2conv_cnn":
        model = GatV2Conv()

    elif model_name == "gineconv_cnn":
        model = GineConv()

    elif model_name == "mpnn_cnn":
        model = MPNNModel()

    elif model_name == "ResGatedGraphConv":
        model = ResGatedModel()
    # elif global_model_name == 'transformer_cnn':
    return model


class CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.input_shape = (26, 1000)

        self.extraction = nn.Sequential(
            nn.Conv1d(in_channels=26, out_channels=64, kernel_size=3),
            nn.ReLU(inplace=True),
            nn.Conv1d(in_channels=64, out_channels=256, kernel_size=3),
            nn.ReLU(inplace=True),
            nn.Conv1d(in_channels=256, out_channels=1024, kernel_size=5),
            nn.ReLU(inplace=True),
            nn.Dropout(0.25),
            nn.AdaptiveMaxPool1d(1),
        )
        self.output = nn.Sequential(
            nn.Linear(self.__get_output_shape(self.input_shape), 768), nn.ReLU()
        )

    def __get_output_shape(self, input_dim):

        x = self.extraction(torch.rand(1, *input_dim))
        # x = x.data.shape
        x = torch.reshape(x, (-1, 1))
        x = x.data.shape[0]

        return x

    def forward(self, input):
        input = input.to(DEVICE, dtype=torch.float)
        x = self.extraction(input)
        # print(f'x.shape before flattening:{x.data.shape}')
        x = x.view(x.size(0), -1)
        # print(f'x.shape:{x.data.shape}')
        x = self.output(x)
        # print(f'x.shape:{x.data.shape}')

        return x


# Drug encodings with transformer
class Transformer(nn.Module):
    def __init__(self):
        super().__init__()
        # Load the two models
        self.drug_model = TransformerConv(in_channels=33, out_channels=256, edge_dim=7)
        self.target_model = CNN()

        # A final multi-layer perception for classification
        self.mlp = nn.Sequential(
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 1),
            nn.ReLU(),
        )

    def forward(self, drug, target):
        # Input the drug and target model
        _drug = self.drug_model(
            x=drug.x, edge_index=drug.edge_index.long(), edge_attr=drug.edge_attr
        )
        _drug = scatter_mean(_drug, drug.batch, dim=0)
        _target = self.target_model(target)

        # Concatenate two features
        _output = torch.cat((_drug, _target), 1)

        output = self.mlp(_output)

        return output


# Drug encodings with gatv2conv
class GatV2Conv(nn.Module):
    def __init__(self):
        super().__init__()
        # Load the two models
        self.drug_model = GATv2Conv(in_channels=33, out_channels=256, edge_dim=7)
        self.target_model = CNN()

        # A final multi-layer perception for classification
        self.mlp = nn.Sequential(
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 1),
            nn.ReLU(),
        )

    def forward(self, drug, target):
        # Input the drug and target model
        _drug = self.drug_model(
            x=drug.x, edge_index=drug.edge_index.long(), edge_attr=drug.edge_attr
        )
        _drug = scatter_mean(_drug, drug.batch, dim=0)
        _target = self.target_model(target)

        # Concatenate two features
        _output = torch.cat((_drug, _target), 1)

        output = self.mlp(_output)

        return output


class GINEEncoder(nn.Module):
    def __init__(self, in_channels, out_channels, e_dim):
        super().__init__()

        self.encode_node = nn.Sequential(
            nn.Linear(in_channels, out_channels), nn.ReLU()
        )

        self.linear1 = nn.Sequential(
            nn.Linear(out_channels, out_channels),
            nn.ReLU(),
            nn.Linear(out_channels, out_channels),
            nn.ReLU(),
        )

        self.linear2 = nn.Sequential(
            nn.Linear(out_channels, out_channels),
            nn.ReLU(),
            nn.Linear(out_channels, out_channels),
            nn.ReLU(),
        )

        self.conv1 = GINEConv(self.linear1, edge_dim=e_dim, train_eps=True)
        self.conv2 = GINEConv(self.linear2, edge_dim=e_dim, train_eps=True)

        self.dropout = nn.Dropout(0.2)
        self.activation = F.relu

    def forward(self, x, edge_index, edge_attr):

        x = self.encode_node(x)
        x = self.dropout(self.activation(self.conv1(x, edge_index, edge_attr)))
        drug_embedding = self.activation(self.conv2(x, edge_index, edge_attr))

        return drug_embedding


# Drug encodings with GINEConv
class GineConv(nn.Module):
    def __init__(self):
        super().__init__()
        # Load the two models
        self.drug_model = GINEEncoder(in_channels=33, out_channels=256, e_dim=7)
        self.target_model = CNN()

        # A final multi-layer perception for classification
        self.mlp = nn.Sequential(
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 1),
            nn.ReLU(),
        )

    def forward(self, drug, target):
        # Input the drug and target model
        _drug = self.drug_model(
            x=drug.x, edge_index=drug.edge_index.long(), edge_attr=drug.edge_attr
        )
        _drug = scatter_mean(_drug, drug.batch, dim=0)
        _target = self.target_model(target)

        # Concatenate two features
        _output = torch.cat((_drug, _target), 1)

        output = self.mlp(_output)

        return output


class MPNN(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super().__init__(aggr="add")  # "Add" aggregation (Step 5).
        self.lin = Linear(in_channels, out_channels, bias=False)
        self.bias = Parameter(torch.Tensor(out_channels))

        self.reset_parameters()

    def reset_parameters(self):
        self.lin.reset_parameters()
        self.bias.data.zero_()

    def forward(self, data):
        x = data.x
        edge_index = data.edge_index

        edge_index = edge_index.long()

        x = self.lin(x)
        # print(f'Linear_Transform_Shape{x.data.shape}')

        # Step 3: Compute normalization.
        row, col = edge_index
        deg = degree(col, x.size(0), dtype=x.dtype)
        deg_inv_sqrt = deg.pow(-0.5)
        deg_inv_sqrt[deg_inv_sqrt == float("inf")] = 0
        norm = deg_inv_sqrt[row] * deg_inv_sqrt[col]

        # Step 4-5: Start propagating messages.
        out = self.propagate(edge_index, x=x, norm=norm)

        # Step 6: Apply a final bias vector.
        out += self.bias

        # print(f'Out_Shape{out.data.shape}')

        out = scatter_mean(out, data.batch, dim=0)

        return out

    def message(self, x_j, norm):
        return norm.view(-1, 1) * x_j


# Drug encodings with mpnn
class MPNNModel(nn.Module):
    def __init__(self):
        super().__init__()
        # Load the two models
        self.drug_model = MPNN(33, 256)
        self.target_model = CNN()

        # A final multi-layer perception for classification
        self.mlp = nn.Sequential(
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 1),
            # nn.ReLU()
        )

    def forward(self, drug, target):
        # Input the drug and target model
        _drug = self.drug_model(drug)
        _target = self.target_model(target)

        # Concatenate two features
        _output = torch.cat((_drug, _target), 1)

        output = self.mlp(_output)

        return output


# Drug encodings with ResGatedGraphConv
class ResGatedModel(nn.Module):
    def __init__(self):
        super().__init__()
        # Load the two models
        self.drug_model = ResGatedGraphConv(in_channels=33, out_channels=256)
        self.target_model = CNN()

        # A final multi-layer perception for classification
        self.mlp = nn.Sequential(
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(1024, 1024),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 1),
            nn.ReLU(),
        )

    def forward(self, drug, target):
        # Input the drug and target model
        _drug = self.drug_model(x=drug.x, edge_index=drug.edge_index.long())
        _drug = scatter_mean(_drug, drug.batch, dim=0)
        _target = self.target_model(target)

        # Concatenate two features
        _output = torch.cat((_drug, _target), 1)

        output = self.mlp(_output)

        return output
