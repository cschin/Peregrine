from setuptools import setup
import os
os.environ["peregrine_base"] = os.path.abspath(os.path.pardir)

setup(name='peregrine_pypy',
      version='0.1',
      install_requires=["networkx==2.4"],
      scripts = ["scripts/ovlp_to_graph.py", "scripts/graph_to_path.py"])
