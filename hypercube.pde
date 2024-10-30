import ddf.minim.*;
import ddf.minim.analysis.*;
import peasy.*;

Minim minim;
AudioOutput out;
FFT fftLeft, fftRight;
PeasyCam cam;
Hypercube cube;

final int SAMPLE_RATE = 44100;
final int BUFFER_SIZE = 1024;

// Sweep parameters
float currentFreq = 50;
float minFreq = 50;
float maxFreq = 1000;
float sweepDuration = 1.0; // seconds
float phaseLeft = 0.0;
float phaseRight = 0.0;
float sweepTime = 0.0;
boolean sweepingUp = true;
float stereoWidth = 0.0; // 0.0 = mono, 1.0 = full stereo

// Time domain parameters
final int TIME_HISTORY = 1000; // milliseconds of history
float lastUpdateTime = 0;
ArrayList<TimeSlice> timeHistory;

void setup() {
    size(800, 800, P3D);
    background(0);
    
    cube = new Hypercube(256, 32);
    timeHistory = new ArrayList<TimeSlice>();
    
    minim = new Minim(this);
    out = minim.getLineOut(Minim.STEREO, BUFFER_SIZE);
    
    // Initialize FFTs for both channels
    fftLeft = new FFT(BUFFER_SIZE, SAMPLE_RATE);
    fftRight = new FFT(BUFFER_SIZE, SAMPLE_RATE);
    fftLeft.window(FFT.HAMMING);
    fftRight.window(FFT.HAMMING);
    
    cube.InitSmoothSpectrum(fftLeft.specSize());

    frameRate(60);

    cam = new PeasyCam(this, 800);
    cam.setMinimumDistance(200);
    cam.setMaximumDistance(2000);
    sphereDetail(6);
    
    // Initialize time history with empty slices
    for (int i = 0; i < cube.NumLayers; i++) {
        timeHistory.add(new TimeSlice(cube.GridSize));
    }
}

void updateSweep() {
    float timeIncrement = 1.0 / frameRate;
    sweepTime += timeIncrement;
    
    if (sweepTime >= sweepDuration) {
        sweepTime = 0;
        sweepingUp = !sweepingUp;
    }
    
    float sweepProgress = sweepTime / sweepDuration;
    if (sweepingUp) {
        currentFreq = minFreq * pow(maxFreq/minFreq, sweepProgress);
        stereoWidth = sweepProgress;
    } else {
        currentFreq = maxFreq * pow(minFreq/maxFreq, sweepProgress);
        stereoWidth = 1.0 - sweepProgress;
    }
}

void generateStereoSweepBuffer(float[] bufferLeft, float[] bufferRight) {
    float phaseIncrement = TWO_PI * currentFreq / SAMPLE_RATE;
    
    for (int i = 0; i < bufferLeft.length; i++) {
        float signal = 0.5 * sin(phaseLeft);
        
        bufferLeft[i] = signal * (1.0 - stereoWidth * 0.5);
        bufferRight[i] = signal * (1.0 + stereoWidth * 0.5);
        
        phaseLeft += phaseIncrement;
        if (phaseLeft > TWO_PI) phaseLeft -= TWO_PI;
    }
}

void updateTimeHistory() {
    float currentTime = millis();
    float deltaTime = currentTime - lastUpdateTime;
    
    if (deltaTime >= (TIME_HISTORY / cube.NumLayers)) {
        // Remove oldest time slice and add new one at the beginning
        timeHistory.remove(timeHistory.size() - 1);
        timeHistory.add(0, new TimeSlice(cube.GridSize));
        
        lastUpdateTime = currentTime;
    }
}

void draw() {
    background(0);
    stroke(255, 0, 255);
    noFill();
    box(cube.CubeSize);


    // Framerate display
    fill(255);
    textSize(16);
    textAlign(RIGHT, TOP);
    text("FPS: " + int(frameRate), cube.CubeSize / 2, cube.CubeSize / 2, cube.CubeSize / 2);

    updateSweep();
    updateTimeHistory();
    
    float[] audioBufferLeft = new float[BUFFER_SIZE];
    float[] audioBufferRight = new float[BUFFER_SIZE];
    generateStereoSweepBuffer(audioBufferLeft, audioBufferRight);
    
    // Process FFTs
    fftLeft.forward(audioBufferLeft);
    fftRight.forward(audioBufferRight);

    translate(-(cube.CubeSize / 2), -(cube.CubeSize / 2), 0);
    noLights();

    // Store current state in newest time slice
    TimeSlice currentSlice = timeHistory.get(0);
    
    for (int layer = 0; layer < cube.NumLayers; layer++) {
        float z = map(layer, 0, cube.NumLayers - 1, -(cube.CubeSize / 2), cube.CubeSize / 2);
        TimeSlice slice = timeHistory.get(layer);

        for (int row = 0; row < cube.GridSize; row++) {
            float rowFreq = cube.YMin * pow(cube.YMax/cube.YMin, (float)row/(cube.GridSize-1));
            float y = map(row, 0, cube.GridSize - 1, cube.CubeSize, 0);
            
            // Calculate frequency proximity factor
            float freqDistance = abs(rowFreq - currentFreq);
            float relativeDistance = freqDistance / currentFreq;
            float freqFactor = (relativeDistance < 0.1) ? map(relativeDistance, 0, 0.1, 1.0, 0.0) : 0.0;

            for (int col = 0; col < cube.GridSize; col++) {
                float x = map(col, 0, cube.GridSize - 1, 0, cube.CubeSize);
                
                float baseSize = 1.0;
                float maxSize = 5.0;
                
                // Calculate stereo position factor (0 = center, 1 = edges)
                float stereoPosition = (float)col / (cube.GridSize - 1);
                float centerDistance = abs(stereoPosition - 0.5) * 2.0; // 0 at center, 1 at edges
                
                // Enhanced stereo factor calculation
                float stereoFactor;
                if (stereoWidth < 0.1) {
                    // Mostly mono: emphasize center
                    stereoFactor = (1.0 - centerDistance) * map(stereoWidth, 0, 0.1, 1.0, 0.8);
                } else if (stereoWidth > 0.9) {
                    // Mostly stereo: emphasize sides
                    stereoFactor = centerDistance * map(stereoWidth, 0.9, 1.0, 0.8, 1.0);
                } else {
                    // Transitional: gradual shift from center to sides
                    float monoComponent = (1.0 - centerDistance) * (1.0 - stereoWidth);
                    float stereoComponent = centerDistance * stereoWidth;
                    stereoFactor = max(monoComponent, stereoComponent);
                }
                
                // Apply frequency factor to determine final particle size
                float particleSize = baseSize;
                if (freqFactor > 0) {
                    particleSize += (maxSize - baseSize) * freqFactor * stereoFactor;
                }
                
                // Store current values in the newest time slice
                if (layer == 0) {
                    currentSlice.setSizeAt(row, col, particleSize);
                    currentSlice.setStereoFactorAt(row, col, stereoFactor);
                }
                
                // Get historical values for rendering
                float historicalSize = slice.getSizeAt(row, col);
                float historicalStereoFactor = slice.getStereoFactorAt(row, col);

                // Color calculation with enhanced contrast
                if (historicalSize > baseSize) {
                    float intensity = map(historicalStereoFactor, 0, 1, 0.2, 1.0);
                    float blue = 255 * intensity;
                    float alpha = 255 * intensity;
                    fill(0, blue, blue, alpha);
                } else {
                    fill(0, 0, 0, 0);
                }
                
                noStroke();
                pushMatrix();
                translate(x, y, z);
                sphere(historicalSize);
                popMatrix();
            }
        }
    }
}

// Class to store a single time slice of the visualization
class TimeSlice {
    private float[][] sizes;
    private float[][] stereoFactors;
    
    public TimeSlice(int gridSize) {
        sizes = new float[gridSize][gridSize];
        stereoFactors = new float[gridSize][gridSize];
        
        // Initialize with default values
        for (int i = 0; i < gridSize; i++) {
            for (int j = 0; j < gridSize; j++) {
                sizes[i][j] = 1.0; // baseSize
                stereoFactors[i][j] = 0.0;
            }
        }
    }
    
    public void setSizeAt(int row, int col, float size) {
        sizes[row][col] = size;
    }
    
    public float getSizeAt(int row, int col) {
        return sizes[row][col];
    }
    
    public void setStereoFactorAt(int row, int col, float factor) {
        stereoFactors[row][col] = factor;
    }
    
    public float getStereoFactorAt(int row, int col) {
        return stereoFactors[row][col];
    }
}

public class Hypercube {
    public int CubeSize = 400;
    public int NumParticles = 256;
    public int GridSize = (int) sqrt(NumParticles);
    public int NumLayers = 21;

    public float YMin = 50.0;    // Hz - Match sweep minimum
    public float YMax = 1000.0;  // Hz - Match sweep maximum
    public float[] smoothedSpectrum;
    
    //public Hyperercube() {
    //}

    public Hypercube(int numParticles, int numLayers) {
      NumParticles = numParticles;
      GridSize = (int) sqrt(numParticles);
      NumLayers = numLayers;
    }

    public void InitSmoothSpectrum(int specSize) {
        smoothedSpectrum = new float[specSize];
    }
}

public static class FFTUtil {
    private static final float MIN_AMPLITUDE = 1e-6f;
    private static final float REF_AMPLITUDE = 1.0f;
    
    public static float toDB(float amplitude) {
        float absAmplitude = Math.abs(amplitude);
        float safeAmplitude = Math.max(absAmplitude, MIN_AMPLITUDE);
        return 20 * (float)Math.log10(safeAmplitude / REF_AMPLITUDE);
    }
}
