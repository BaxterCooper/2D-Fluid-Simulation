#include <math.h>
#include <SFML/Graphics.hpp>

int N = 128;

int SCALE = 4;
int ITERATIONS = 16;

float TIMESTEP = 0.1f;
float DIFFUSION = 0;
float VISCOSITY = 0.0000001f;

#define IX(x, y) ((x) + (y) * N)

struct FluidSquare {
    int size;
    float dt;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;
    float *Vz;

    float *Vx0;
    float *Vy0;
    float *Vz0;
};

FluidSquare *FluidSquareCreate(int size, float diffusion, float viscosity, float dt)
{
    FluidSquare *square = (FluidSquare *) malloc(sizeof(*square));
    int N = size;
    
    square->size = size;
    square->dt = dt;
    square->diff = diffusion;
    square->visc = viscosity;
    
    square->s = (float *) calloc(N * N * N, sizeof(float));
    square->density = (float *) calloc(N * N * N, sizeof(float));
    
    square->Vx = (float *) calloc(N * N * N, sizeof(float));
    square->Vy = (float *) calloc(N * N * N, sizeof(float));
    
    square->Vx0 = (float *) calloc(N * N * N, sizeof(float));
    square->Vy0 = (float *) calloc(N * N * N, sizeof(float));
    
    return square;
}

static void setBounds(int b, float *x)
{
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)]   = b == 2 ? -x[IX(i, 1)]   : x[IX(i, 1)];
        x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
    }
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)]   = b == 1 ? -x[IX(1, j)]   : x[IX(1, j)];
        x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
    }

    x[IX(0, 0)]     = 0.5 * (x[IX(1, 0)]     + x[IX(0, 1)]);
    x[IX(0, N-1)]   = 0.5 * (x[IX(1, N-1)]   + x[IX(0, N-2)]);
    x[IX(N-1, 0)]   = 0.5 * (x[IX(N-2, 0)]   + x[IX(N-1, 1)]);
    x[IX(N-1, N-1)] = 0.5 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);
}

static void linearSolve(int b, float *x, float *x0, float a, float c)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < ITERATIONS; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)]
                        + a*(    x[IX(i+1, j)]
                                +x[IX(i-1, j)]
                                +x[IX(i, j+1)]
                                +x[IX(i, j-1)]
                        )) * cRecip;
            }
        }
        setBounds(b, x);
    }
}

static void diffuse (int b, float *x, float *x0, float diff, float dt)
{
    float a = dt * diff * (N - 2) * (N - 2);
    linearSolve(b, x, x0, a, 1 + 6 * a);
}

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float dt)
{
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    
    for(int j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
        for(int i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x    = ifloat - tmp1; 
            y    = jfloat - tmp2;
            
            if(x < 0.5f) x = 0.5f; 
            if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
            i0 = floorf(x); 
            i1 = i0 + 1.0f;
            if(y < 0.5f) y = 0.5f; 
            if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
            j0 = floorf(y);
            j1 = j0 + 1.0f;
            
            s1 = x - i0; 
            s0 = 1.0f - s1; 
            t1 = y - j0; 
            t0 = 1.0f - t1;
            
            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;
            
            d[IX(i, j)] =
                s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
        }
    }

    setBounds(b, d);
}

static void project(float *velocX, float *velocY, float *p, float *div)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f*(
                        velocX[IX(i+1, j)]
                    -velocX[IX(i-1, j)]
                    +velocY[IX(i  , j+1)]
                    -velocY[IX(i  , j-1)]
                )/N;
            p[IX(i, j)] = 0;
        }
    }
    setBounds(0, div); 
    setBounds(0, p);
    linearSolve(0, p, div, 1, 6);
    
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
        }
    }

    setBounds(1, velocX);
    setBounds(2, velocY);
}

void FluidSquareStep(FluidSquare *square)
{
    int N          = square->size;
    float visc     = square->visc;
    float diff     = square->diff;
    float dt       = square->dt;
    float *Vx      = square->Vx;
    float *Vy      = square->Vy;
    float *Vz      = square->Vz;
    float *Vx0     = square->Vx0;
    float *Vy0     = square->Vy0;
    float *Vz0     = square->Vz0;
    float *s       = square->s;
    float *density = square->density;
    
    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);
    
    project(Vx0, Vy0, Vx, Vy);
    
    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);
    
    project(Vx, Vy, Vx0, Vy0);
    
    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
}

void FluidSquareAddDensity(FluidSquare *square, int x, int y, float amount)
{
    square->density[IX(x, y)] += amount;
}

void FluidSquareAddVelocity(FluidSquare *square, int x, int y, float amountX, float amountY)
{
    square->Vx[IX(x, y)] += amountX;
    square->Vy[IX(x, y)] += amountY;
}

void FluidSquareFadeDensity(FluidSquare *square) {
    for (int i = 0; i < N * N; i++) {
        if (square->density[i] > 0) {
            square->density[i] -= 0.01;
        }
    }
}

void render(FluidSquare *square, sf::RenderWindow &window) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int x = i * SCALE;
            int y = j * SCALE;
            float density = square->density[IX(i, j)];

            sf::RectangleShape cell;
            cell.setSize(sf::Vector2f(SCALE, SCALE));
            cell.setPosition(i * SCALE, j * SCALE);
            cell.setFillColor(sf::Color(density, density, density));

            window.draw(cell);
        }
    }
}

int main()
{
    FluidSquare *square = FluidSquareCreate(N, DIFFUSION, VISCOSITY, TIMESTEP);

    sf::RenderWindow window(sf::VideoMode(N*SCALE, N*SCALE), "2D Fluid Simulation");

    sf::Vector2i previousMouse = sf::Mouse::getPosition(window);
	sf::Vector2i currentMouse = sf::Mouse::getPosition(window);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::Black);

        currentMouse = sf::Mouse::getPosition(window);

        FluidSquareAddDensity(square, currentMouse.x / SCALE, currentMouse.y / SCALE, 200);
        FluidSquareAddVelocity(square, currentMouse.x / SCALE, currentMouse.y / SCALE, currentMouse.x - previousMouse.x, currentMouse.y - previousMouse.y);
        FluidSquareFadeDensity(square);

        render(square, window);

        FluidSquareStep(square);

        previousMouse = currentMouse;
        
        window.display();
    }

    return 0;
}
