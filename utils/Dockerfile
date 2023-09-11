# Use an official Python runtime as a base image
FROM python:3.8

# Set the working directory in the container to /app
WORKDIR /findadapt

# Copy the current directory contents into the container at /app
ADD . /findadapt

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir pyahocorasick

# Make port 80 available to the world outside this container
EXPOSE 8080

# Define environment variable

# Set the default command for the container
CMD ["-h"]

# Specify the entrypoint of your application
ENTRYPOINT ["/findadapt/findadapt"]
