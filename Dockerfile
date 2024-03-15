# Use an official Python runtime as a parent image
FROM python:3.8

# Set the working directory inside the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Clone the OpenPNM repository and install it
RUN git clone https://github.com/PMEAL/OpenPNM \
    && cd OpenPNM \
    && git checkout v3.0.0 \
    && pip install --no-cache-dir -r requirements.txt \
    && pip install -e .
    
# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements_PNM.txt

# RUN cp patch/* OpenPNM/ -r
# Patch required when running pressuredrop fitting scripts only

# Run the main.py file with parameters
CMD ["python", "PNM_main_V3_Final_OpenPNM.py"]
